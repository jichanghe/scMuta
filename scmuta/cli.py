"""scmuta call: BAM(s) -> candidate sites -> per-cell calls with calibration.

Plate data (one BAM per cell, e.g. Smart-seq2):
  scmuta call --plate cells.tsv --fasta GRCh38.fa --regions exons.bed --negative ctrl_cells.txt --out out/
Droplet data (one BAM, cell barcode tag):
  scmuta call --bam possorted.bam --tag CB --cells barcodes.tsv --fasta GRCh38.fa --regions exons.bed \\
              --stranded --out out/
Known sites only (e.g. an engineered mutation), skipping discovery:
  scmuta call --plate cells.tsv --fasta GRCh38.fa --sites known.tsv --truth labels.csv --out out/
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np
import pandas as pd

from . import __version__
from .model import SiteFilters
from .pileup import ReadFilter, count, discover, plate_bams
from .pipeline import run_counts


def _regions(path):
    bed = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], comment="#")
    return [(str(c), int(s), int(e)) for c, s, e in bed.itertuples(index=False)]


def _lines(path):
    with open(path) as fh:
        return [x.strip().split("\t")[0].split(",")[0] for x in fh if x.strip() and not x.startswith("#")]


def main(argv=None):
    ap = argparse.ArgumentParser(prog="scmuta", description="Calibrated per-cell SNV calling")
    ap.add_argument("--version", action="version", version=f"scMuta {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("call")
    src = c.add_mutually_exclusive_group(required=True)
    src.add_argument("--plate", help="table: cell<TAB>bam (one BAM per cell)")
    src.add_argument("--bam", help="one BAM with a cell-barcode tag")
    c.add_argument("--tag", default="CB")
    c.add_argument("--cells", help="file with the cell barcodes to keep (droplet)")
    c.add_argument("--fasta", required=True)
    grp = c.add_mutually_exclusive_group(required=True)
    grp.add_argument("--regions", help="BED of regions to scan (e.g. expressed exons)")
    grp.add_argument("--sites", help="known sites: chrom, pos (0-based), ref, alt; skips discovery")
    c.add_argument("--negative", help="file listing negative-control cells (cannot carry the variants)")
    c.add_argument("--target-fp", type=float, default=0.001, help="max false-call rate in negative controls")
    c.add_argument("--threshold", type=float, default=None, help="fixed posterior threshold (overrides calibration)")
    c.add_argument("--min-alt-reads", type=int, default=3)
    c.add_argument("--min-alt-cells", type=int, default=2)
    c.add_argument("--min-mapq", type=int, default=20)
    c.add_argument("--min-baseq", type=int, default=20)
    c.add_argument("--umi-tag", default="UB", help="UMI tag for --dedup umi (droplet)")
    c.add_argument("--dedup", choices=["none", "umi", "start"], default=None,
                   help="count molecules, not reads (default: umi for --bam, none for --plate)")
    c.add_argument("--stranded", action="store_true", help="stranded protocol (10x): no strand-bias filter")
    c.add_argument("--editing-bed", help="known RNA-editing sites (chrom, pos)")
    c.add_argument("--flag-all-ag", action="store_true", help="flag every A>G/T>C as possible editing")
    c.add_argument("--truth", help="optional cell,label CSV (label 1 = carries the variant) for known sites")
    c.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    rf = ReadFilter(min_mapq=a.min_mapq, min_baseq=a.min_baseq)
    if a.plate:
        bams, tag, keep = plate_bams(a.plate), None, None
    else:
        bams, tag = a.bam, a.tag
        keep = set(_lines(a.cells)) if a.cells else None
    if a.sites:
        sites = pd.read_csv(a.sites, sep=None, engine="python")
        sites.columns = [x.lower() for x in sites.columns]
        sites = sites[["chrom", "pos", "ref", "alt"]]
    else:
        sites = discover(bams, a.fasta, _regions(a.regions), rf, tag, keep, a.min_alt_reads, a.min_alt_cells)
    if sites.empty:
        raise SystemExit("no candidate sites found")
    dedup = a.dedup or ("umi" if a.bam else "none")
    A, D, s, cells = count(bams, a.fasta, sites, rf, tag, keep, umi_tag=a.umi_tag, dedup=dedup)
    neg = None
    if a.negative:
        negset = set(_lines(a.negative))
        neg = np.array([c in negset for c in cells])
    filters = SiteFilters(stranded_protocol=a.stranded, editing_bed=a.editing_bed, flag_all_ag=a.flag_all_ag)
    res = run_counts(A, D, s, cells, negative=neg, filters=filters, target_fp=a.target_fp, threshold=a.threshold,
                     min_alt_cells=1 if a.sites else a.min_alt_cells, out=a.out)
    if a.truth:
        tr = pd.read_csv(a.truth).set_index(pd.read_csv(a.truth).columns[0]).iloc[:, 0]
        called = set(res["calls"].cell)
        lab = tr.reindex(cells)
        cov = (D.sum(1) > 0)
        df = pd.DataFrame({"cell": cells, "label": lab.to_numpy(), "called": [c in called for c in cells], "covered": cov})
        df.to_csv(os.path.join(a.out, "truth_check.tsv"), sep="\t", index=False)
        d = df[df.covered & df.label.notna()]
        m = {"carriers_covered": int((d.label == 1).sum()), "carriers_called": int(((d.label == 1) & d.called).sum()),
             "noncarriers_covered": int((d.label == 0).sum()), "noncarriers_called": int(((d.label == 0) & d.called).sum())}
        print(json.dumps(m))
    print(json.dumps(res["summary"], default=float, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
