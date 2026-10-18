"""Pilot for 10x single-nucleus RNA from several donors (e.g. GSE174332, STARsolo output).

Answers two questions before a full project:
  1. coverage: per cell type, how many genome positions per nucleus have >= 1 / >= 2 molecules, and how many
     positions are covered (>= 2 molecules) in >= k nuclei of one type in one donor (what a clone needs)
  2. scMuta: candidate SNVs on the chosen chromosomes, each donor tested against the other donors' nuclei
     as negative controls (a somatic mutation of donor A cannot be in donor B); germline variants are judged
     within each donor, because a heterozygous SNP of donor A is absent from donor B and would otherwise
     look somatic

Usage:
  python examples/snrna_pilot.py --donor 214:STARsolo/214 --donor 301:STARsolo/301 \
      --fasta GRCh38.fa --chroms chr22 --out pilot_out
Each donor folder is a STARsolo output (Aligned.sortedByCoord.out.bam(.bai), Solo.out/GeneFull/filtered/).
"""
import argparse
import gzip
import json
import os
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam
from scipy import io, sparse

from scmuta.model import SiteFilters
from scmuta.pileup import ReadFilter, count, discover
from scmuta.pipeline import run_counts_chunked
from scipy import sparse as sp
from scmuta.util import log

MARKERS = {
    "ExN": ["SLC17A7", "SATB2", "RBFOX3", "SNAP25"],
    "InN": ["GAD1", "GAD2", "SST", "PVALB", "VIP"],
    "Oligo": ["MBP", "PLP1", "MOG", "ST18"],
    "OPC": ["PDGFRA", "VCAN", "OLIG1"],
    "Astro": ["AQP4", "GFAP", "SLC1A2", "ALDH1L1"],
    "Micro": ["P2RY12", "CSF1R", "CX3CR1", "C3", "DOCK8"],
    "Endo": ["CLDN5", "FLT1", "VWF"],
}


def _open(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_solo(folder):
    """Barcodes, genes and matrix (cells x genes) from STARsolo GeneFull/filtered."""
    base = os.path.join(folder, "Solo.out", "GeneFull", "filtered")
    f = lambda n: next(os.path.join(base, x) for x in (n, n + ".gz") if os.path.exists(os.path.join(base, x)))
    bc = [l.strip() for l in _open(f("barcodes.tsv"))]
    genes = [l.rstrip("\n").split("\t")[1] for l in _open(f("features.tsv"))]
    m = io.mmread(f("matrix.mtx")).tocsr().T.tocsr()     # genes x cells -> cells x genes
    return bc, genes, m


def cell_types(m, genes):
    """Marker-score label per nucleus (pilot quality; a full study would cluster)."""
    lib = np.asarray(m.sum(1)).ravel()
    norm = sparse.diags(1e4 / np.maximum(lib, 1)) @ m
    norm.data = np.log1p(norm.data)
    gi = {g: i for i, g in enumerate(genes)}
    scores = {}
    for t, gs in MARKERS.items():
        idx = [gi[g] for g in gs if g in gi]
        scores[t] = np.asarray(norm[:, idx].mean(1)).ravel() if idx else np.zeros(m.shape[0])
    sc = pd.DataFrame(scores)
    lab = sc.idxmax(1).to_numpy()
    lab[sc.max(1).to_numpy() < 0.2] = "unassigned"
    return lab, lib


def coverage(bam, cells, chroms, step=10, chunk=1_000_000, rf=None):
    """Per cell: positions (every `step`-th base, scaled back) with >=1 and >=2 UMIs; per position: cells
    with >=2 UMIs. Streams the sorted BAM chunk by chunk."""
    rf = rf or ReadFilter()
    ci = {c: i for i, c in enumerate(cells)}
    ge1 = np.zeros(len(cells))
    ge2 = np.zeros(len(cells))
    pos_cells2 = []           # per sampled position: array of cell indices with >= 2 UMIs
    with pysam.AlignmentFile(bam) as f:
        for chrom in chroms:
            L = f.get_reference_length(chrom)
            for c0 in range(0, L, chunk):
                c1 = min(L, c0 + chunk)
                umis = defaultdict(set)      # (pos, cell) -> UMIs
                for r in f.fetch(chrom, c0, c1):
                    if not rf.read_ok(r) or not r.has_tag("CB") or not r.has_tag("UB"):
                        continue
                    j = ci.get(r.get_tag("CB"))
                    if j is None:
                        continue
                    u = r.get_tag("UB")
                    for b0, b1 in r.get_blocks():
                        lo = max(b0, c0)
                        lo += (-lo) % step
                        for p in range(lo, min(b1, c1), step):
                            umis[(p, j)].add(u)
                bypos = defaultdict(list)
                for (p, j), us in umis.items():
                    ge1[j] += 1
                    if len(us) >= 2:
                        ge2[j] += 1
                        bypos[p].append(j)
                pos_cells2.extend(np.array(v) for v in bypos.values())
            log(f"coverage {os.path.basename(os.path.dirname(bam))} {chrom} done")
    return ge1 * step, ge2 * step, pos_cells2


def main():
    import functools, builtins
    global print
    print = functools.partial(builtins.print, flush=True)   # results appear even if the run is interrupted
    ap = argparse.ArgumentParser()
    ap.add_argument("--donor", action="append", required=True, help="NAME:STARsolo_folder (repeat)")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--chroms", default="chr22")
    ap.add_argument("--min-cells-type", type=int, default=30)
    ap.add_argument("--out", default="snrna_pilot")
    ap.add_argument("--editing-bed", default=None, help="known RNA-editing sites (e.g. REDIportal hg38), chrom pos")
    ap.add_argument("--min-wt-nuclei", type=int, default=5,
                    help="a somatic call needs >= this many covered (>=2 molecules) wild-type nuclei of the same donor")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    chroms = a.chroms.split(",")
    donors = [x.split(":", 1) for x in a.donor]

    # 1. cells and cell types
    meta, bams = [], {}
    for d, folder in donors:
        bc, genes, m = read_solo(folder)
        lab, lib = cell_types(m, genes)
        meta.append(pd.DataFrame({"donor": d, "barcode": bc, "cell": [f"{d}:{b}" for b in bc],
                                  "cell_type": lab, "umis": lib}))
        bams[d] = os.path.join(folder, "Aligned.sortedByCoord.out.bam")
    meta = pd.concat(meta, ignore_index=True)
    print("\nNuclei per donor and cell type:")
    print(pd.crosstab(meta.cell_type, meta.donor, margins=True).to_string())
    print("\nMedian UMIs per nucleus:")
    print(meta.pivot_table(index="cell_type", columns="donor", values="umis", aggfunc="median").round(0).to_string())

    # 2. coverage
    cov_rows, clone_rows = [], []
    for d, _ in donors:
        md = meta[meta.donor == d].reset_index(drop=True)
        g1, g2, pc = coverage(bams[d], md.barcode.tolist(), chroms)
        md["pos_ge1"], md["pos_ge2"] = g1, g2
        cov_rows.append(md)
        types = md.cell_type.to_numpy()
        for t in sorted(set(types)):
            n_t = int((types == t).sum())
            if n_t < a.min_cells_type:
                continue
            counts = np.array([(types[v] == t).sum() for v in pc]) if pc else np.zeros(0)
            clone_rows.append({"donor": d, "cell_type": t, "nuclei": n_t,
                               **{f"positions_ge2umi_in>={k}_nuclei": int((counts >= k).sum()) * 10
                                  for k in (2, 5, 10, 20)}})
    cov = pd.concat(cov_rows, ignore_index=True)
    cov.to_csv(os.path.join(a.out, "coverage_per_nucleus.tsv"), sep="\t", index=False)
    clone = pd.DataFrame(clone_rows)
    clone.to_csv(os.path.join(a.out, "clone_detectable_positions.tsv"), sep="\t", index=False)
    print(f"\nPositions per nucleus on {a.chroms} with >=1 / >=2 molecules (median):")
    print(cov.groupby(["donor", "cell_type"])[["pos_ge1", "pos_ge2"]].median().round(0).to_string())
    print("\nPositions covered (>=2 molecules) in at least k nuclei of one type (what a clone needs):")
    print(clone.to_string(index=False))

    # 3. scMuta: discover per donor, union, count per donor, stack
    rf = ReadFilter()
    regions = []
    with pysam.FastaFile(a.fasta) as fa:
        for c in chroms:
            L = fa.get_reference_length(c)
            regions += [(c, s, min(L, s + 5_000_000)) for s in range(0, L, 5_000_000)]
    cands = []
    for d, _ in donors:
        keep = set(meta.barcode[meta.donor == d])
        cd = discover(bams[d], a.fasta, regions, rf, "CB", keep, min_alt_reads=3, min_alt_cells=2)
        cands.append(cd[["chrom", "pos", "ref", "alt"]])
    sites = pd.concat(cands).drop_duplicates().sort_values(["chrom", "pos"]).reset_index(drop=True)
    log(f"{len(sites)} candidate alleles (union over donors)")
    As, Ds, cells_all, feats = [], [], [], None
    for d, _ in donors:
        keep = set(meta.barcode[meta.donor == d])
        A, D, s, cells = count(bams[d], a.fasta, sites, rf, "CB", keep, umi_tag="UB", dedup="umi")
        # count() returns only cells with reads at candidate sites, in sorted order; align to all nuclei (sparse:
        # each nucleus covers few sites, so dense cells x sites matrices would not fit in memory)
        A, D = sp.csr_matrix(A), sp.csr_matrix(D)
        full = meta.barcode[meta.donor == d].tolist()
        idx = {c: i for i, c in enumerate(cells)}
        rows = np.array([idx.get(c, -1) for c in full])
        P = sp.csr_matrix((np.ones((rows >= 0).sum()), (np.flatnonzero(rows >= 0), rows[rows >= 0])),
                          shape=(len(full), len(cells)))
        As.append((P @ A).tocsr())
        Ds.append((P @ D).tocsr())
        del A, D
        cells_all += [f"{d}:{c}" for c in full]
        num = ["alt_fwd", "alt_rev", "ref_fwd", "ref_rev", "alt_near_end", "alt_near_junction", "alt_near_clip",
               "alt_baseq_sum"]
        if feats is None:
            feats = s.copy()
        else:
            feats[num] = feats[num].to_numpy() + s[num].to_numpy()
    A = sp.vstack(As).tocsr()
    D = sp.vstack(Ds).tocsr()
    log(f"counts: {A.shape[0]} nuclei x {A.shape[1]} sites, {D.nnz} covered pairs ({D.nnz / np.prod(D.shape):.2e} filled)")
    cm = meta.set_index("cell").loc[cells_all]
    donor_of = cm.donor.to_numpy()

    # germline within each donor. A germline SNP is in every cell type; a somatic clone is normally confined to
    # one lineage (e.g. microglia), even when large. So: alt in >= 30% of covered nuclei (>= 3 covered) in at
    # least 2 cell types, or, when only one type is covered, alt in >= 50% of covered nuclei with VAF >= 0.3.
    types_of = cm.cell_type.to_numpy()
    germ = {}
    for d, _ in donors:
        m = donor_of == d
        n_types_alt = np.zeros(len(sites), int)
        n_types_cov = np.zeros(len(sites), int)
        for t in set(types_of[m]) - {"unassigned"}:
            mt = m & (types_of == t)
            cov_t = np.asarray((D[mt] > 0).sum(0)).ravel()
            alt_t = np.asarray((A[mt] > 0).sum(0)).ravel()
            ok = cov_t >= 3
            n_types_cov += ok
            n_types_alt += ok & (alt_t >= 0.3 * np.maximum(cov_t, 1))
        a_, d_ = np.asarray(A[m].sum(0)).ravel(), np.asarray(D[m].sum(0)).ravel()
        cov_, alt_ = np.asarray((D[m] > 0).sum(0)).ravel(), np.asarray((A[m] > 0).sum(0)).ravel()
        single = (n_types_cov <= 1) & (alt_ >= 0.5 * np.maximum(cov_, 1)) & (a_ >= 0.3 * np.maximum(d_, 1)) & (cov_ >= 3)
        germ[d] = (n_types_alt >= 2) | single
        log(f"donor {d}: {germ[d].sum()} germline-like alleles (heterozygous/homozygous SNPs)")

    results, summary = [], {"candidates": int(len(sites))}
    for d, _ in donors:
        keep_sites = ~germ[d]
        neg = donor_of != d
        res = run_counts_chunked(A[:, keep_sites], D[:, keep_sites], feats[keep_sites].reset_index(drop=True), cells_all,
                         negative=neg, filters=SiteFilters(stranded_protocol=True, editing_cluster_window=50, editing_bed=a.editing_bed), target_fp=0.001,
                         min_alt_cells=2, out=os.path.join(a.out, f"scmuta_{d}"))
        st = res["sites"]
        calls = res["calls"].merge(cm[["cell_type"]], left_on="cell", right_index=True)
        calls["donor"] = d
        # informativeness: a somatic variant must be absent from covered nuclei of the same donor; otherwise a
        # barely covered germline heterozygous SNP looks identical
        m = (donor_of == d)
        Dk, Ak = D[:, keep_sites][m], A[:, keep_sites][m]
        cov2 = (Dk >= 2)
        wt = np.asarray(cov2.sum(0)).ravel() - np.asarray(cov2.multiply(Ak > 0).sum(0)).ravel()
        st["donor_wt_nuclei"] = wt
        st["somatic_supported"] = st["variant"] & (st["donor_wt_nuclei"] >= a.min_wt_nuclei)
        ok_sites = set(zip(st.chrom[st.somatic_supported], st.pos[st.somatic_supported]))
        calls["somatic_supported"] = [(c, p) in ok_sites for c, p in zip(calls.chrom, calls.pos)]
        results.append(calls)
        # burden per genotypable cell-site pair (>= 2 molecules at PASS sites), by cell type: coverage-normalised
        pass_idx = np.flatnonzero((st["filter"] == "PASS").to_numpy())
        Dp = D[:, keep_sites][:, pass_idx]
        gen_pairs = np.asarray((Dp[m] >= 2).sum(1)).ravel()
        tm = types_of[m]
        bur = []
        for t in sorted(set(tm)):
            n_pairs = int(gen_pairs[tm == t].sum())
            n_calls = int(((calls.cell_type == t) & calls.somatic_supported).sum())
            bur.append({"cell_type": t, "nuclei": int((tm == t).sum()), "genotypable_pairs": n_pairs,
                        "supported_calls": n_calls, "calls_per_1e5_pairs": round(1e5 * n_calls / max(n_pairs, 1), 2),
                        "all_calls": int((calls.cell_type == t).sum())})
        bur = pd.DataFrame(bur)
        bur.to_csv(os.path.join(a.out, f"burden_{d}.tsv"), sep="\t", index=False)
        summary[d] = {"germline_like": int(germ[d].sum()), **{k: res["summary"][k] for k in
                      ("pass_filters", "variants", "cell_calls", "threshold")},
                      "somatic_supported_variants": int(st["somatic_supported"].sum()),
                      "removed_by_negative_controls": int((st["filter"] == "in_negative_controls").sum()),
                      "filters": res["summary"]["filters"]}
        print(f"\n=== donor {d} vs other donors as negatives")
        print(json.dumps(summary[d], indent=1, default=float))
        if res["summary"]["variants"]:
            v = st[st.variant][["chrom", "pos", "ref", "alt", "subst", "vaf_pooled", "cells_called", "donor_wt_nuclei",
                                "somatic_supported", "negative_cells_called"]]
            print(v.assign(pos=v.pos + 1).head(40).to_string(index=False))
        print("burden by cell type (coverage-normalised):")
        print(bur.to_string(index=False))
    pd.concat(results).to_csv(os.path.join(a.out, "somatic_calls.tsv"), sep="\t", index=False)
    json.dump(summary, open(os.path.join(a.out, "summary.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
