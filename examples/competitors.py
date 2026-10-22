"""Head-to-head on the cell-line benchmark: SComatic, scAllele, scMuta and common-practice rules, scored against the
same truth (examples/cellline_benchmark.py must have run first; it writes the truth/test split and counts).

Every tool sees only the TEST-half cells. Each run records wall time, CPU time and peak memory.

  python examples/competitors.py prep     --bench bench_cellline --star STAR
  python examples/competitors.py scomatic --bench bench_cellline --fasta GRCh38.fa --scomatic-dir tools/SComatic
  python examples/competitors.py scallele --bench bench_cellline --fasta GRCh38.fa --star STAR
  python examples/competitors.py scmuta   --bench bench_cellline --fasta GRCh38.fa --star STAR
  python examples/competitors.py score    --bench bench_cellline --fasta GRCh38.fa --star STAR

Tools run inside their own conda environments through `conda run -n <env>` (--env to change the name).
Per-cell calls:
  SComatic : PASS sites from BaseCellCalling.step2; a cell is called when it belongs to a cell type the site was
             called in and has >= 1 alt read there (what SingleCellGenotype.py reports with --alt_flag Alt).
             SComatic needs cell types, so it receives the cell-line labels: an advantage the other tools lack.
  scAllele : run per cell (one BAM per cell, as for Smart-seq); a cell is called at SNVs in its VCF with FILTER
             PASS (or '.').
  scMuta   : bench/scmuta/cell_calls.tsv from cellline_benchmark.py; the `scmuta` subcommand reruns the same
             analysis on the test cells only, to measure time and memory.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import resource
import subprocess
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pysam
from scipy import sparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cellline_benchmark as cb  # noqa: E402

from scmuta.model import SiteFilters  # noqa: E402
from scmuta.pileup import ReadFilter, count, discover  # noqa: E402
from scmuta.pipeline import run_counts_chunked  # noqa: E402

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


# --------------------------------------------------------------------------- helpers
def load_split(bench, star):
    t = pd.read_csv(os.path.join(bench, "cells_split.tsv"), sep="\t", dtype=str)
    bams = {}
    for r in t.run_accession:
        hit = glob.glob(os.path.join(star, r, "*.sortedByCoord.out.bam"))
        if hit:
            bams[r] = hit[0]
    return t, bams


def test_cells(t):
    return t[(t.half == "test") & t["cell line"].isin(cb.CASE_LINES + cb.CONTROL_LINES)]


class Timer:
    """Wall, CPU and peak RSS of this process and its children."""

    def __enter__(self):
        self.t0 = time.time()
        self.c0 = resource.getrusage(resource.RUSAGE_CHILDREN)
        self.s0 = resource.getrusage(resource.RUSAGE_SELF)
        return self

    def __exit__(self, *exc):
        c = resource.getrusage(resource.RUSAGE_CHILDREN)
        s = resource.getrusage(resource.RUSAGE_SELF)
        self.wall = time.time() - self.t0
        self.cpu = (c.ru_utime + c.ru_stime - self.c0.ru_utime - self.c0.ru_stime +
                    s.ru_utime + s.ru_stime - self.s0.ru_utime - self.s0.ru_stime)
        self.peak_gb = max(c.ru_maxrss, s.ru_maxrss) / 1024 / 1024     # ru_maxrss is in KB on Linux


def save_resources(out, tm, extra=None):
    os.makedirs(out, exist_ok=True)
    d = {"wall_hours": tm.wall / 3600, "cpu_hours": tm.cpu / 3600, "peak_memory_gb": tm.peak_gb}
    d.update(extra or {})
    json.dump(d, open(os.path.join(out, "resources.json"), "w"), indent=1)
    print(json.dumps(d, indent=1), flush=True)


def sh(cmd, log):
    print("+", cmd if isinstance(cmd, str) else " ".join(cmd), flush=True)
    with open(log, "a") as fh:
        subprocess.run(cmd, shell=isinstance(cmd, str), check=True, stdout=fh, stderr=subprocess.STDOUT)


# --------------------------------------------------------------------------- prep: one CB-tagged BAM for SComatic
def prep(a):
    t, bams = load_split(a.bench, a.star)
    tc = test_cells(t)
    out = os.path.join(a.bench, "competitors", "prep")
    tag_dir = os.path.join(out, "tagged")
    os.makedirs(tag_dir, exist_ok=True)
    parts = []
    tags_seen = set()
    for c in tc.run_accession:
        dst = os.path.join(tag_dir, f"{c}.bam")
        parts.append(dst)
        if os.path.exists(dst):
            continue
        with pysam.AlignmentFile(bams[c]) as src, pysam.AlignmentFile(dst + ".tmp", "wb", template=src) as o:
            for k, r in enumerate(src.fetch(until_eof=True)):
                if k < 1000:
                    tags_seen |= {x for x, _ in r.get_tags()}
                r.set_tag("CB", c)
                o.write(r)
        os.replace(dst + ".tmp", dst)
    if tags_seen and not {"nM", "NH"} <= tags_seen:
        print(f"WARNING: SComatic filters on nM and NH tags; BAM tags seen: {sorted(tags_seen)}")
    merged = os.path.join(out, "test_cells.bam")
    pysam.merge("-f", "-@", str(a.threads), merged, *parts)
    pysam.index(merged)
    pd.DataFrame({"Index": tc.run_accession, "Cell_type": tc["cell line"]}).to_csv(
        os.path.join(out, "meta.tsv"), sep="\t", index=False)
    print(f"merged BAM: {merged} ({len(parts)} cells); meta.tsv written", flush=True)


# --------------------------------------------------------------------------- SComatic
def scomatic(a):
    d = os.path.abspath(a.scomatic_dir)
    prepd = os.path.abspath(os.path.join(a.bench, "competitors", "prep"))
    out = os.path.abspath(os.path.join(a.bench, "competitors", "scomatic"))
    for sub in ("split", "counts", "calling", "tmp"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    log = os.path.join(out, "scomatic.log")
    py = f"conda run --no-capture-output -n {a.env} python"
    fa = os.path.abspath(a.fasta)
    with Timer() as tm:
        sh(f"cd {d} && {py} scripts/SplitBam/SplitBamCellTypes.py --bam {prepd}/test_cells.bam "
           f"--meta {prepd}/meta.tsv --id tseq --n_trim 5 --max_nM 5 --max_NH 1 --outdir {out}/split", log)
        for bam in sorted(glob.glob(f"{out}/split/tseq.*.bam")):
            sh(f"cd {d} && {py} scripts/BaseCellCounter/BaseCellCounter.py --bam {bam} --ref {fa} --chrom all "
               f"--out_folder {out}/counts --min_bq 30 --tmp_dir {out}/tmp --nprocs {a.threads}", log)
        sh(f"cd {d} && {py} scripts/MergeCounts/MergeBaseCellCounts.py --tsv_folder {out}/counts "
           f"--outfile {out}/calling/tseq.BaseCellCounts.AllCellTypes.tsv", log)
        sh(f"cd {d} && {py} scripts/BaseCellCalling/BaseCellCalling.step1.py "
           f"--infile {out}/calling/tseq.BaseCellCounts.AllCellTypes.tsv --outfile {out}/calling/tseq --ref {fa}", log)
        sh(f"cd {d} && {py} scripts/BaseCellCalling/BaseCellCalling.step2.py "
           f"--infile {out}/calling/tseq.calling.step1.tsv --outfile {out}/calling/tseq "
           f"--editing {d}/RNAediting/AllEditingSites.hg38.txt --pon {d}/PoNs/PoN.scRNAseq.hg38.tsv", log)
    save_resources(out, tm, {"tool": "SComatic", "threads": a.threads})


def read_scomatic(bench, fasta):
    f = os.path.join(bench, "competitors", "scomatic", "calling", "tseq.calling.step2.tsv")
    if not os.path.exists(f):
        return None
    # SComatic writes "##" metadata lines first, then a "#CHROM ..." header
    with open(f) as fh:
        skip = 0
        for line in fh:
            if line.startswith("##"):
                skip += 1
            else:
                break
    df = pd.read_csv(f, sep="\t", dtype=str, skiprows=skip)
    df.columns = [c.lstrip("#") for c in df.columns]
    df = df[df.FILTER == "PASS"]
    fa = pysam.FastaFile(fasta)
    rows = []
    for _, r in df.iterrows():
        start = int(r.Start)
        # coordinate convention checked against the reference base
        pos0 = start - 1 if fa.fetch(r.CHROM, start - 1, start).upper() == r.REF.upper() else start
        for alt in str(r.ALT).split(","):
            if len(alt) == 1 and len(r.REF) == 1:
                for ct in str(r.Cell_types).split(","):
                    rows.append((r.CHROM, pos0, r.REF.upper(), alt.upper(), ct.strip()))
    return pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt", "cell_type"])


# --------------------------------------------------------------------------- scAllele (per cell)
def scallele(a):
    t, bams = load_split(a.bench, a.star)
    tc = test_cells(t)
    out = os.path.abspath(os.path.join(a.bench, "competitors", "scallele"))
    os.makedirs(out, exist_ok=True)
    log = os.path.join(out, "scallele.log")
    fa = os.path.abspath(a.fasta)
    regions = [x for x in (a.region_chroms or "").split(",") if x]
    with Timer() as tm:
        for c in tc.run_accession:
            if os.path.exists(os.path.join(out, f"{c}.vcf")):
                continue
            if not regions:
                sh(f"conda run --no-capture-output -n {a.env} scAllele -b {os.path.abspath(bams[c])} -g {fa} "
                   f"-o {out}/{c} -n {a.threads}", log)
                continue
            # restricted run (scAllele is slow genome-wide): one call per chromosome, VCFs merged per cell
            parts = []
            for ch in regions:
                pre = os.path.join(out, "parts", f"{c}.{ch}")
                os.makedirs(os.path.dirname(pre), exist_ok=True)
                if not os.path.exists(pre + ".vcf"):
                    sh(f"conda run --no-capture-output -n {a.env} scAllele -b {os.path.abspath(bams[c])} -g {fa} "
                       f"-o {pre} -n {a.threads} -c {ch}", log)
                parts.append(pre + ".vcf")
            with open(os.path.join(out, f"{c}.vcf.tmp"), "w") as o:
                for k, pth in enumerate(parts):
                    for line in open(pth):
                        if line.startswith("#") and k > 0:
                            continue
                        o.write(line)
            os.replace(os.path.join(out, f"{c}.vcf.tmp"), os.path.join(out, f"{c}.vcf"))
    save_resources(out, tm, {"tool": "scAllele", "threads": a.threads, "cells": len(tc),
                             "chromosomes": regions or "all"})


def read_scallele(bench, cells):
    out = os.path.join(bench, "competitors", "scallele")
    rows = []
    found = 0
    for c in cells:
        f = os.path.join(out, f"{c}.vcf")
        if not os.path.exists(f):
            continue
        found += 1
        for line in open(f):
            if line.startswith("#"):
                continue
            x = line.rstrip("\n").split("\t")
            chrom, pos, ref, alts, filt = x[0], int(x[1]), x[3].upper(), x[4].upper(), x[6]
            if filt not in ("PASS", "."):
                continue
            for alt in alts.split(","):
                if len(ref) == 1 and len(alt) == 1:
                    rows.append((c, chrom, pos - 1, ref, alt))
    if not found:
        return None
    return pd.DataFrame(rows, columns=["cell", "chrom", "pos", "ref", "alt"])


# --------------------------------------------------------------------------- scMuta timed on test cells
_G = {}


def _disc(args):
    chrom, length = args
    return discover(_G["bams"], _G["fasta"], [(chrom, 0, length)], ReadFilter(), None, None,
                    min_alt_reads=3, min_alt_cells=2)


def _cnt(sites):
    return count(_G["bams"], _G["fasta"], sites, ReadFilter(), None, None, dedup="start")


def scmuta_timed(a):
    t, bams = load_split(a.bench, a.star)
    tc = test_cells(t)
    cells = tc.run_accession.tolist()
    _G.update(bams={c: bams[c] for c in cells}, fasta=a.fasta)
    fa = pysam.FastaFile(a.fasta)
    chroms = [c for c in a.chroms.split(",") if c in fa.references]
    out = os.path.join(a.bench, "competitors", "scmuta_timed")
    with Timer() as tm:
        with Pool(a.threads, maxtasksperchild=1) as p:
            cand = pd.concat(p.map(_disc, [(c, fa.get_reference_length(c)) for c in chroms]))
        cand = cand.sort_values(["chrom", "pos", "alt"]).reset_index(drop=True)
        tasks = []
        for c in chroms:
            g = cand[cand.chrom == c].reset_index(drop=True)
            tasks += [g.iloc[b:b + a.chunk_sites].reset_index(drop=True) for b in range(0, len(g), a.chunk_sites)]
        with Pool(a.threads, maxtasksperchild=1) as p:
            parts = list(p.imap(_cnt, tasks))
        As, Ds, ss = [], [], []
        for A, D, s, cl in parts:
            ci = {c: k for k, c in enumerate(cl)}
            A2 = np.zeros((len(cells), A.shape[1]), np.int32)
            D2 = np.zeros_like(A2)
            for k, c in enumerate(cells):
                if c in ci:
                    A2[k], D2[k] = A[ci[c]], D[ci[c]]
            As.append(sparse.csc_matrix(A2))
            Ds.append(sparse.csc_matrix(D2))
            ss.append(s)
        A, D = sparse.hstack(As).tocsc(), sparse.hstack(Ds).tocsc()
        s = pd.concat(ss).reset_index(drop=True)
        neg = tc["cell line"].isin(cb.CONTROL_LINES).to_numpy()
        r = run_counts_chunked(A, D, s, cells, negative=neg, filters=SiteFilters(germline_cell_frac=0.95),
                       target_fp=a.target_fp, min_alt_cells=2, out=out)
    save_resources(out, tm, {"tool": "scMuta", "threads": a.threads, "cells": len(cells),
                             "variants": r["summary"]["variants"]})


# --------------------------------------------------------------------------- scoring
def _cnt_all(args):
    sites = args
    A, D, s, cl = count(_G["bams"], _G["fasta"], sites, ReadFilter(), None, None, dedup="start")
    order = _G["cells"]
    ci = {c: k for k, c in enumerate(cl)}
    A2 = np.zeros((len(order), len(sites)), np.int32)
    D2 = np.zeros_like(A2)
    for k, c in enumerate(order):
        if c in ci:
            A2[k], D2[k] = A[ci[c]], D[ci[c]]
    return sparse.csr_matrix(A2), sparse.csr_matrix(D2), s


def score(a):
    t, bams = load_split(a.bench, a.star)
    cells = t.run_accession.tolist()
    z = np.load(os.path.join(a.bench, "counts.npz"), allow_pickle=True)
    A, D = sparse.csr_matrix(z["A"].item()), sparse.csr_matrix(z["D"].item())
    sites = pd.read_csv(os.path.join(a.bench, "sites_counted.tsv"), sep="\t")
    line = t["cell line"].to_numpy()
    truth = (t.half == "truth").to_numpy()
    tc_idx = np.flatnonzero((~truth) & np.isin(line, cb.CASE_LINES + cb.CONTROL_LINES))
    tcells = [cells[i] for i in tc_idx]

    calls = {}
    sm = pd.read_csv(os.path.join(a.bench, "scmuta", "cell_calls.tsv"), sep="\t")
    calls["scMuta"] = sm[["cell", "chrom", "pos", "alt"]]
    sc = read_scomatic(a.bench, a.fasta) if "SComatic" not in (a.skip or "") else None
    if sc is not None:
        lt = dict(zip(t.run_accession, t["cell line"]))
        calls["SComatic"] = sc   # expanded to cells below (needs counts)
    sa = read_scallele(a.bench, tcells) if "scAllele" not in (a.skip or "") else None
    if sa is not None:
        calls["scAllele"] = sa[["cell", "chrom", "pos", "alt"]]

    # count any called site not yet in the matrix, on all cells (truth half too, for its label)
    key = set(zip(sites.chrom, sites.pos, sites.alt))
    extra = []
    for name, df in calls.items():
        for c_, p_, a_ in zip(df.chrom, df.pos, df.alt):
            if (c_, int(p_), a_) not in key:
                extra.append((c_, int(p_), a_))
                key.add((c_, int(p_), a_))
    if extra:
        fa = pysam.FastaFile(a.fasta)
        ex = pd.DataFrame(extra, columns=["chrom", "pos", "alt"])
        ex["ref"] = [fa.fetch(c, p, p + 1).upper() for c, p in zip(ex.chrom, ex.pos)]
        ex = ex[["chrom", "pos", "ref", "alt"]].sort_values(["chrom", "pos"]).reset_index(drop=True)
        print(f"counting {len(ex)} sites called by competitors but not in the scMuta candidate set", flush=True)
        _G.update(bams={c: bams[c] for c in cells}, fasta=a.fasta, cells=cells)
        with Pool(a.threads, maxtasksperchild=1) as p:
            res = list(p.imap(_cnt_all, [g.iloc[b:b + a.chunk_sites].reset_index(drop=True)
                                         for c in ex.chrom.unique() for g in [ex[ex.chrom == c].reset_index(drop=True)]
                                         for b in range(0, len(g), a.chunk_sites)]))
        A = sparse.hstack([A] + [r[0] for r in res]).tocsr()
        D = sparse.hstack([D] + [r[1] for r in res]).tocsr()
        sites = pd.concat([sites] + [r[2].assign(found_in="competitor") for r in res]).reset_index(drop=True)

    lab, _ = cb.truth_labels(A, D, line, truth)
    if a.chroms_scored:
        keep_ch = set(a.chroms_scored.split(","))
        lab = np.where(sites.chrom.isin(keep_ch).to_numpy(), lab, "excluded")
        print(f"scoring restricted to {sorted(keep_ch)}: {int((lab != 'excluded').sum())} of {len(lab)} sites")
    At, Dt = A[tc_idx].toarray(), D[tc_idx].toarray()
    lt_arr = line[tc_idx]
    case = np.isin(lt_arr, cb.CASE_LINES)
    neg = np.isin(lt_arr, cb.CONTROL_LINES)
    sidx = {k: j for j, k in enumerate(zip(sites.chrom, sites.pos.astype(int), sites.alt))}
    cidx = {c: i for i, c in enumerate(tcells)}

    def matrix(df):
        M = np.zeros(At.shape, bool)
        for c_, ch, p_, a_ in zip(df.cell, df.chrom, df.pos, df.alt):
            i, j = cidx.get(c_), sidx.get((ch, int(p_), a_))
            if i is not None and j is not None:
                M[i, j] = True
        return M

    results = []
    for name, df in calls.items():
        if name == "SComatic":
            M = np.zeros(At.shape, bool)
            for ch, p_, a_, ct in zip(df.chrom, df.pos, df.alt, df.cell_type):
                j = sidx.get((ch, int(p_), a_))
                if j is not None:
                    M[:, j] |= (lt_arr == ct) & (At[:, j] >= 1)
        else:
            M = matrix(df)
        results.append(cb.metrics(M, At, Dt, lt_arr, lab, case, name))
    tested = sites.found_in.fillna("").str.contains("test").to_numpy()
    vaf = np.where(Dt > 0, At / np.maximum(Dt, 1), 0)
    n2 = (At >= 2) & (vaf >= 0.1) & tested[None, :]
    n2 &= (n2 & case[:, None]).sum(0)[None, :] >= 2
    results.append(cb.metrics(n2, At, Dt, lt_arr, lab, case, ">=2 alt, VAF>=0.1"))
    results.append(cb.metrics(n2 & ~((At[neg] >= 1).any(0))[None, :], At, Dt, lt_arr, lab, case,
                              ">=2 alt, VAF>=0.1, drop sites in controls"))
    tab = pd.DataFrame(results).set_index("method")
    for tool, sub in (("scMuta", "scmuta_timed"), ("SComatic", "scomatic"), ("scAllele", "scallele")):
        f = os.path.join(a.bench, "competitors", sub, "resources.json")
        if os.path.exists(f) and tool in tab.index:
            r = json.load(open(f))
            for k in ("wall_hours", "cpu_hours", "peak_memory_gb"):
                tab.loc[tool, k] = r[k]
    out = os.path.join(a.bench, "competitors")
    tag = ("_" + a.chroms_scored.replace(",", "_")) if a.chroms_scored else ""
    os.makedirs(out, exist_ok=True)
    tab.to_csv(os.path.join(out, f"head_to_head{tag}.tsv"), sep="\t")
    pd.set_option("display.width", 220)
    print(tab.T.to_string())
    print("\nnote: SComatic received the cell-line labels as cell types (its model needs them); the other tools "
          "did not. Pairs are test-half cells; truth from the held-out truth half.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["prep", "scomatic", "scallele", "scmuta", "score"])
    ap.add_argument("--bench", default="bench_cellline")
    ap.add_argument("--star", default="STAR")
    ap.add_argument("--fasta")
    ap.add_argument("--scomatic-dir")
    ap.add_argument("--env")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--chroms", default=",".join(CHROMS))
    ap.add_argument("--chunk-sites", type=int, default=20000)
    ap.add_argument("--region-chroms", default=None, help="scallele: run only on these chromosomes (comma list)")
    ap.add_argument("--chroms-scored", default=None, help="score: only sites on these chromosomes count")
    ap.add_argument("--skip", default=None, help="score: tools to leave out, e.g. scAllele")
    ap.add_argument("--target-fp", type=float, default=0.001)
    a = ap.parse_args()
    if a.env is None:
        a.env = {"scomatic": "SComatic", "scallele": "scallele"}.get(a.step, "")
    {"prep": prep, "scomatic": scomatic, "scallele": scallele, "scmuta": scmuta_timed, "score": score}[a.step](a)


if __name__ == "__main__":
    main()
