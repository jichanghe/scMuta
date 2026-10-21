"""Genome-wide per-cell benchmark on plate scRNA-seq of mixed cell lines (TARGET-seq GSE105451 validation set).

Ground truth comes from cell identity, not from any caller:
  * cells of each (cell line, method) stratum are split at random into a TRUTH half and a TEST half
  * truth labels use only the TRUTH half, pooled per line:
      private to line L : pooled VAF in L >= --truth-vaf with >= 6 alt reads in >= 2 cells, and essentially
                          absent (alt <= max(1, 0.5% of depth), depth >= 20) in every other line.
                          These are L's germline heterozygous SNPs and L's somatic mutations (e.g. JAK2 V617F in
                          SET2): every L cell carries them, no other cell does.
      absent            : depth >= 20 and essentially no alt reads in every line -> any call is false.
      in_controls       : also carried by the control individual (HSPC): shared germline SNP, RNA editing or a
                          recurrent artefact -> not a somatic variant of the case cells; calls there are counted.
      everything else   : ambiguous, not scored.
  * callers see only the TEST half: the case cells are SET2 + JURKAT mixed together (a site private to one line
    behaves like a clonal somatic variant present in a subset of cells), the HSPC cells are the negative
    controls (a different individual, as in a cohort).
Callers compared: scMuta, and common practice rules (>= 2 alt reads and VAF >= 0.1; the same plus dropping sites
seen in the controls; any alt read). Metrics per cell x site pair: sensitivity at carrier pairs, false-positive
rate at non-carrier pairs and at absent sites, precision, and for scMuta the confident wild-type calls and the
controls-estimated false-positive rate against the observed one.

  python examples/cellline_benchmark.py --cells metadata/cells.tsv --star STAR --fasta GRCh38.fa \
      --threads 8 --out bench_cellline
"""
from __future__ import annotations

import argparse
import glob
import json
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pysam
from scipy import sparse

from scmuta.model import SiteFilters, genotype_calls
from scmuta.pileup import ReadFilter, count, discover
from scmuta.pipeline import run_counts_chunked

CASE_LINES = ("SET2", "JURKAT")
CONTROL_LINES = ("HSPC",)
_G = {}


def _discover(args):
    chrom, length, which = args
    bams = _G["bams_" + which]
    d = discover(bams, _G["fasta"], [(chrom, 0, length)], ReadFilter(), None, None,
                 min_alt_reads=3, min_alt_cells=2)
    return d


def _count(args):
    chrom, sites = args
    A, D, s, cells = count(_G["bams_all"], _G["fasta"], sites, ReadFilter(), None, None, dedup="start")
    order = _G["cells"]
    ci = {c: k for k, c in enumerate(cells)}
    A2 = np.zeros((len(order), len(sites)), np.int32)
    D2 = np.zeros_like(A2)
    for k, c in enumerate(order):
        if c in ci:
            A2[k], D2[k] = A[ci[c]], D[ci[c]]
    return sparse.csr_matrix(A2), sparse.csr_matrix(D2), s


def truth_labels(A, D, line, truth, vaf=0.15, min_dep=20, absent_frac=0.005):
    """Per site: 'private_<L>', 'absent' or 'ambiguous' from the truth-half cells only."""
    lines = sorted(set(line))
    stats = {}
    for L in lines:
        m = truth & (line == L)
        a = np.asarray(A[m].sum(0)).ravel()
        d = np.asarray(D[m].sum(0)).ravel()
        cells_alt = np.asarray((A[m] >= 1).sum(0)).ravel()
        stats[L] = (a, d, cells_alt)
    n = A.shape[1]
    lab = np.array(["ambiguous"] * n, dtype=object)
    clean = {L: (stats[L][1] >= min_dep) & (stats[L][0] <= np.maximum(1, absent_frac * stats[L][1])) for L in lines}
    absent = np.ones(n, bool)
    for L in lines:
        absent &= clean[L]
    lab[absent] = "absent"
    for L in lines:
        a, d, ca = stats[L]
        present = (d > 0) & (a >= 6) & (ca >= 2) & (a >= vaf * np.maximum(d, 1))
        others = np.ones(n, bool)
        for M in lines:
            if M != L:
                others &= clean[M]
        lab[present & others] = f"private_{L}"
    # present in the control individual too: shared germline SNP, RNA editing or a recurrent artefact; whatever it
    # is, it is not a somatic variant of the case cells, so a somatic call there is a false somatic call
    for C in CONTROL_LINES:
        if C in stats:
            a, d, ca = stats[C]
            lab[(lab == "ambiguous") & (a >= 3) & (ca >= 2) & (a >= 0.05 * np.maximum(d, 1))] = "in_controls"
    return lab, stats


def metrics(called, A, D, line, lab, case, name, wt=None):
    """Pair-level metrics on case cells; carrier = cell's line owns the private site."""
    rows = {}
    cov = (D > 0)
    priv = np.array([s.startswith("private_") and s[8:] in CASE_LINES for s in lab])
    owner = np.array([s[8:] if s.startswith("private_") else "" for s in lab])
    carrier = (line[:, None] == owner[None, :]) & priv[None, :]
    noncarrier = (~carrier) & priv[None, :]
    absent = (lab == "absent")[None, :]
    cc = case[:, None]
    tp = (called & carrier & cov & cc).sum()
    fn = ((~called) & carrier & cov & cc).sum()
    fp_nc = (called & noncarrier & cov & cc).sum()
    fp_ab = (called & absent & cov & cc).sum()
    n_nc = (noncarrier & cov & cc).sum()
    n_ab = (absent & cov & cc).sum()
    rows["carrier_pairs_covered"] = int(tp + fn)
    rows["sensitivity"] = float(tp / max(tp + fn, 1))
    rows["fpr_noncarrier"] = float(fp_nc / max(n_nc, 1))
    rows["fpr_absent_sites"] = float(fp_ab / max(n_ab, 1))
    rows["true_calls"] = int(tp)
    rows["false_calls"] = int(fp_nc + fp_ab)
    rows["precision"] = float(tp / max(tp + fp_nc + fp_ab, 1))
    rows["false_variant_sites"] = int((called & absent & cc).any(0).sum() + 0)
    ctl = (lab == "in_controls")[None, :]
    rows["somatic_calls_at_sites_present_in_controls"] = int((called & ctl & cc).sum())
    rows["sites_present_in_controls_reported"] = int((called & ctl & cc).any(0).sum())
    rows["calls_at_ambiguous_sites"] = int((called & (lab == "ambiguous")[None, :] & cc).sum())
    # sensitivity by depth of the cell at the site
    for lo, hi in ((1, 2), (3, 5), (6, 10), (11, 10 ** 9)):
        m = carrier & cc & (D >= lo) & (D <= hi)
        rows[f"sens_depth_{lo}-{hi if hi < 10 ** 9 else 'inf'}"] = float((called & m).sum() / max(m.sum(), 1))
    if wt is not None:
        w = wt & cc
        rows["wild_type_calls_correct"] = int((w & noncarrier & cov).sum() + (w & absent & cov).sum())
        rows["wild_type_calls_at_carriers"] = int((w & carrier & cov).sum())
        rows["carrier_pairs_wild_type_rate"] = float((w & carrier & cov).sum() / max((carrier & cov & cc).sum(), 1))
    rows["method"] = name
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", required=True, help="cells.tsv: run_accession, cell line, method")
    ap.add_argument("--star", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--methods", default="SMART-seq+,TARGET-seq", help="library methods to use")
    ap.add_argument("--chroms", default=",".join([f"chr{i}" for i in range(1, 23)] + ["chrX"]))
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--target-fp", type=float, default=0.001)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--chunk-sites", type=int, default=20000,
                    help="sites per counting task (memory per worker grows with this)")
    ap.add_argument("--out", default="bench_cellline")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    t = pd.read_csv(a.cells, sep="\t", dtype=str)
    t = t[t["method"].isin(a.methods.split(","))]
    bams = {}
    for r in t.run_accession:
        hit = glob.glob(os.path.join(a.star, r, "*.sortedByCoord.out.bam"))
        if hit and os.path.exists(hit[0] + ".bai"):
            bams[r] = hit[0]
    t = t[t.run_accession.isin(bams)].reset_index(drop=True)
    rng = np.random.default_rng(a.seed)
    t["half"] = "test"
    for _, g in t.groupby(["cell line", "method"]):
        idx = g.index.to_numpy().copy()
        rng.shuffle(idx)
        t.loc[idx[: len(idx) // 2], "half"] = "truth"
    t.to_csv(os.path.join(a.out, "cells_split.tsv"), sep="\t", index=False)
    print(t.groupby(["cell line", "method", "half"]).size().unstack(fill_value=0), flush=True)

    cells = t.run_accession.tolist()
    _G.update(fasta=a.fasta, cells=cells, bams_all={c: bams[c] for c in cells},
              bams_truth={c: bams[c] for c in t.run_accession[t.half == "truth"]},
              bams_test={c: bams[c] for c in t.run_accession[t.half == "test"]})
    fa = pysam.FastaFile(a.fasta)
    chroms = [c for c in a.chroms.split(",") if c in fa.references]
    lens = {c: fa.get_reference_length(c) for c in chroms}

    cache = os.path.join(a.out, "counts.npz")
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        A, D = sparse.csr_matrix(z["A"].item()), sparse.csr_matrix(z["D"].item())
        sites = pd.read_csv(os.path.join(a.out, "sites_counted.tsv"), sep="\t")
        print(f"loaded cached counts: {A.shape}", flush=True)
    else:
        cand_file = os.path.join(a.out, "candidates.tsv")
        if os.path.exists(cand_file):
            found = pd.read_csv(cand_file, sep="\t")
            print(f"loaded cached candidates: {len(found)}", flush=True)
        else:
            with Pool(a.threads, maxtasksperchild=1) as pool:
                disc = pool.map(_discover, [(c, lens[c], w) for c in chroms for w in ("truth", "test")])
            tags = [w for c in chroms for w in ("truth", "test")]
            cand = pd.concat([d.assign(found_in=w) for d, w in zip(disc, tags)])
            found = cand.groupby(["chrom", "pos", "ref", "alt"]).found_in.agg(
                lambda x: ",".join(sorted(set(x)))).reset_index()
            found = found.sort_values(["chrom", "pos", "alt"]).reset_index(drop=True)
            found.to_csv(cand_file, sep="\t", index=False)
        print(f"candidates: {len(found)} ({(found.found_in.str.contains('test')).sum()} found in test cells)", flush=True)
        # count in small blocks of sites: the per-worker tallies grow with sites x cells
        tasks = []
        for c in chroms:
            g = found[found.chrom == c].reset_index(drop=True)
            for b0 in range(0, len(g), a.chunk_sites):
                tasks.append((c, g.iloc[b0:b0 + a.chunk_sites].reset_index(drop=True)))
        print(f"counting in {len(tasks)} blocks of <= {a.chunk_sites} sites", flush=True)
        with Pool(a.threads, maxtasksperchild=1) as pool:
            res = list(pool.imap(_count, tasks))
        A = sparse.hstack([r[0] for r in res]).tocsr()
        D = sparse.hstack([r[1] for r in res]).tocsr()
        sites = pd.concat([r[2] for r in res]).reset_index(drop=True)
        if "found_in" not in sites.columns:
            sites = sites.merge(found[["chrom", "pos", "alt", "found_in"]], on=["chrom", "pos", "alt"], how="left")
        np.savez(cache, A=np.array(A, dtype=object), D=np.array(D, dtype=object))
        sites.to_csv(os.path.join(a.out, "sites_counted.tsv"), sep="\t", index=False)

    line = t["cell line"].to_numpy()
    truth = (t.half == "truth").to_numpy()
    test = ~truth
    lab, _ = truth_labels(A, D, line, truth)
    sites["truth"] = lab
    print("truth labels (sites):", pd.Series(lab).value_counts().to_dict(), flush=True)

    # ---- callers on TEST cells only, at sites discovered in test cells
    tsites = np.flatnonzero(sites.found_in.fillna("").str.contains("test").to_numpy())
    tc = np.flatnonzero(test & np.isin(line, CASE_LINES + CONTROL_LINES))
    At_sp = A[tc][:, tsites].tocsc()
    Dt_sp = D[tc][:, tsites].tocsc()
    At, Dt = At_sp.toarray(), Dt_sp.toarray()
    lt = line[tc]
    labt = lab[tsites]
    case = np.isin(lt, CASE_LINES)
    neg = np.isin(lt, CONTROL_LINES)
    st = sites.iloc[tsites][["chrom", "pos", "ref", "alt"] + [c for c in sites.columns if c.startswith(("alt_", "ref_"))
                                                            or c in ("context", "repeat_extension")]]
    st = st.reset_index(drop=True)
    # case cells are a mix of two lines, so a line-private variant is in ~50% of case cells: germline is decided
    # by the controls design here, not by the in-case cell fraction
    filt = SiteFilters(germline_cell_frac=0.95)
    r = run_counts_chunked(At_sp, Dt_sp, st, [cells[i] for i in tc], negative=neg, filters=filt,
                           target_fp=a.target_fp, min_alt_cells=2, out=os.path.join(a.out, "scmuta"))
    thr = r["summary"]["threshold"]
    var = r["sites"]["variant"].to_numpy()
    pcs, ps = r["posterior_call_sparse"], r["posterior_sparse"]
    called_sm = (pcs >= thr).toarray() & var[None, :]
    # confident wild type: the cell's reads favour wild type over the site's mutant model by >= 10 (as genotype_calls)
    post = ps.toarray()
    covered = post > 0
    pi = np.clip(r["sites"]["frac_mutant"].to_numpy(float), 1e-6, 1 - 1e-6)[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        bf_wt = pi * (1 - post) / ((1 - pi) * post)
    wt_sm = covered & (np.nan_to_num(bf_wt, nan=0, posinf=np.inf) >= 10) & ~called_sm
    del post, bf_wt

    vaf = np.where(Dt > 0, At / np.maximum(Dt, 1), 0)
    naive2 = (At >= 2) & (vaf >= 0.1)
    naive2 &= (naive2 & case[:, None]).sum(0)[None, :] >= 2
    naive2_ctrl = naive2 & ~((At[neg] >= 1).any(0))[None, :]
    naive1 = (At >= 1)
    naive1 &= (naive1 & case[:, None]).sum(0)[None, :] >= 2

    res = [metrics(called_sm, At, Dt, lt, labt, case, "scMuta", wt=wt_sm),
           metrics(naive2, At, Dt, lt, labt, case, ">=2 alt, VAF>=0.1"),
           metrics(naive2_ctrl, At, Dt, lt, labt, case, ">=2 alt, VAF>=0.1, drop sites in controls"),
           metrics(naive1, At, Dt, lt, labt, case, ">=1 alt read")]
    # the same calls scored separately for each library protocol (calls are made once, on all test cells together)
    meth = t["method"].to_numpy()[tc]
    by_proto = []
    for mname in sorted(set(meth)):
        cm = case & (meth == mname)
        for nm, M in (("scMuta", called_sm), (">=2 alt, VAF>=0.1", naive2),
                      (">=2 alt, VAF>=0.1, drop sites in controls", naive2_ctrl)):
            row = metrics(M, At, Dt, lt, labt, cm, nm)
            row["protocol"] = mname
            row["case_cells"] = int(cm.sum())
            by_proto.append(row)
    by_proto = pd.DataFrame(by_proto)
    by_proto.to_csv(os.path.join(a.out, "benchmark_by_protocol.tsv"), sep="\t", index=False)
    print("\nby protocol:\n" + by_proto[["protocol", "method", "case_cells", "carrier_pairs_covered", "sensitivity",
                                          "fpr_noncarrier", "precision", "false_calls"]].to_string(index=False),
          flush=True)
    # scMuta across posterior thresholds (precision / sensitivity trade-off at matched sensitivity), and where the
    # control (HSPC) cells' calls fall: if most are at sites the controls' own individual carries, the controls-based
    # error estimate is inflated by their germline variants, not by sequencing / RT errors
    usable = (r["sites"]["filter"] == "PASS").to_numpy()
    curve, negbd = [], []
    for t_ in (0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999):
        ct = (pcs >= t_).toarray()
        v_ = usable & ((ct & case[:, None]).sum(0) >= 2)
        m = metrics(ct & v_[None, :], At, Dt, lt, labt, case, f"scMuta@{t_}")
        curve.append(m)
        nn = ct[neg][:, usable]
        lab_u = labt[usable]
        row = {"threshold": t_, "control_calls": int(nn.sum()),
               "control_pairs": int((Dt[neg][:, usable] > 0).sum())}
        for L in ("absent", "private_HSPC", "in_controls", "ambiguous", "private_SET2", "private_JURKAT"):
            row[f"at_{L}"] = int(nn[:, lab_u == L].sum())
        negbd.append(row)
    curve = pd.DataFrame(curve).set_index("method")
    curve.to_csv(os.path.join(a.out, "scmuta_threshold_curve.tsv"), sep="\t")
    negbd = pd.DataFrame(negbd)
    negbd.to_csv(os.path.join(a.out, "control_calls_by_truth.tsv"), sep="\t", index=False)
    print("\nscMuta by threshold:\n" + curve[["sensitivity", "fpr_noncarrier", "fpr_absent_sites", "precision",
                                              "false_variant_sites"]].to_string(), flush=True)
    print("\ncontrol-cell calls by what the truth half says about the site:\n" + negbd.to_string(index=False), flush=True)
    tab = pd.DataFrame(res).set_index("method")
    tab.to_csv(os.path.join(a.out, "benchmark.tsv"), sep="\t")
    pd.set_option("display.width", 200)
    print("\n" + tab.T.to_string(), flush=True)
    est = r["summary"].get("held_out_negative_fp_rate")
    sm = tab.loc["scMuta"]
    print(f"\nscMuta calibration: false-positive rate estimated from controls {est}, observed in test cells "
          f"{sm.fpr_noncarrier:.2e} (non-carrier) / {sm.fpr_absent_sites:.2e} (absent sites); threshold {thr}")

    # the known driver as a sanity check
    jak = sites.iloc[tsites].reset_index(drop=True)
    k = np.flatnonzero(((jak.chrom == "chr9") & (jak.pos == 5073769) & (jak.alt == "T")).to_numpy())
    if len(k):
        k = k[0]
        pk = ps[:, k].toarray().ravel()
        pk[pk == 0] = np.nan
        pck = pcs[:, k].toarray().ravel()
        pck[np.isnan(pk)] = np.nan
        gk = genotype_calls(pk[:, None], pck[:, None], thr, r["sites"]["frac_mutant"].iloc[k])[:, 0]
        print(f"\nJAK2 V617F: truth label {labt[k]}; scMuta variant {var[k]}; calls by line:",
              pd.crosstab(lt[case], gk[case]).to_dict())
    json.dump({"benchmark": res, "scmuta_summary": r["summary"],
               "truth_sites": pd.Series(lab).value_counts().to_dict(),
               "test_cells": {L: int((lt == L).sum()) for L in sorted(set(lt))}},
              open(os.path.join(a.out, "benchmark.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
