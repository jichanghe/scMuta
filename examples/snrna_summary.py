#!/usr/bin/env python
"""Aggregate the per-chromosome snrna_pilot runs into the numbers and tables for the motor cortex figure.

The result of the scan is that no somatic variant passed in any donor. That is only informative with a
detection limit attached, so this computes what a clone would have had to look like to be seen:

A variant is reported only if two or more nuclei of the same donor carry it with >= 2 molecules. A clone
occupying a fraction f of a cell type is a random subset of those nuclei, so at a position covered in k
nuclei the clone is expected in f*k of them. Detection needs two, i.e. roughly k >= 2/f. The pilot already
counts, per donor and cell type, how many genomic positions are covered at >= 2 molecules in at least k
nuclei (k = 2, 5, 10, 20), which maps to clone fractions of 100%, 40%, 20% and 10%.

  python examples/snrna_summary.py --dir ~/snrna_genome --out snrna_figure4
"""
from __future__ import annotations

import argparse
import glob
import json
import os

import pandas as pd

K_TO_FRACTION = {2: 1.00, 5: 0.40, 10: 0.20, 20: 0.10}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="folder holding one subfolder per chromosome")
    ap.add_argument("--out", default="snrna_figure4")
    a = ap.parse_args()

    chroms = sorted(d for d in os.listdir(a.dir) if os.path.exists(os.path.join(a.dir, d, "summary.json")))
    print(f"{len(chroms)} chromosomes with results: {' '.join(chroms)}")

    # ---- calls and candidates per donor, summed over chromosomes
    rows = []
    for c in chroms:
        s = json.load(open(os.path.join(a.dir, c, "summary.json")))
        for d, v in s.items():
            if not isinstance(v, dict):
                continue
            rows.append({"chrom": c, "donor": d, "candidates": s.get("candidates", 0),
                         "pass_filters": v.get("pass_filters", 0), "variants": v.get("variants", 0),
                         "cell_calls": v.get("cell_calls", 0),
                         "somatic_supported": v.get("somatic_supported_variants", 0),
                         "germline_like": v.get("germline_like", 0),
                         "removed_by_controls": v.get("removed_by_negative_controls", 0),
                         "threshold": v.get("threshold")})
    calls = pd.DataFrame(rows)
    per_donor = calls.groupby("donor")[["pass_filters", "variants", "cell_calls", "somatic_supported",
                                        "germline_like", "removed_by_controls"]].sum()
    per_donor["threshold_range"] = calls.groupby("donor").threshold.agg(lambda x: f"{x.min()}-{x.max()}")
    total_candidates = calls.groupby("chrom").candidates.first().sum()
    print(f"\ncandidate alleles genome-wide (union over donors): {total_candidates:,}")
    print("\nper donor, summed over chromosomes:")
    print(per_donor.to_string())
    per_donor.to_csv(a.out + "_per_donor.tsv", sep="\t")

    # ---- filter tally: what removed the candidates
    tally = {}
    for c in chroms:
        s = json.load(open(os.path.join(a.dir, c, "summary.json")))
        for d, v in s.items():
            if isinstance(v, dict):
                for k, n in (v.get("filters") or {}).items():
                    tally[k] = tally.get(k, 0) + n
    filt = pd.Series(tally).sort_values(ascending=False)
    print("\nsite filters (summed over donors and chromosomes):")
    print(filt.to_string())
    filt.to_csv(a.out + "_filters.tsv", sep="\t", header=["sites"])

    # ---- detection limit: positions where a clone of a given size could have been found
    cl = pd.concat([pd.read_csv(f, sep="\t").assign(chrom=f.split(os.sep)[-2])
                    for f in glob.glob(os.path.join(a.dir, "*", "clone_detectable_positions.tsv"))],
                   ignore_index=True)
    kcols = {int(c.split(">=")[1].split("_")[0]): c for c in cl.columns if c.startswith("positions_ge2umi_in")}
    summed = cl.groupby(["donor", "cell_type"])[list(kcols.values())].sum()
    nuclei = cl.groupby(["donor", "cell_type"]).nuclei.max()
    det = summed.join(nuclei)
    for k, c in kcols.items():
        det[f"Mb_clone_{int(K_TO_FRACTION[k] * 100)}pct"] = (summed[c] / 1e6).round(2)
    det = det[["nuclei"] + [f"Mb_clone_{int(K_TO_FRACTION[k] * 100)}pct" for k in sorted(kcols, reverse=True)]]
    print("\nmegabases where a clone of the given size in that cell type could have been detected:")
    print(det.to_string())
    det.to_csv(a.out + "_detection_limit.tsv", sep="\t")

    # pooled over donors, for the figure
    pooled = det.groupby("cell_type").agg({"nuclei": "sum",
                                           **{c: "mean" for c in det.columns if c.startswith("Mb_clone")}})
    pooled = pooled.sort_values("nuclei", ascending=False).round(2)
    print("\npooled over donors (nuclei summed, megabases averaged per donor):")
    print(pooled.to_string())
    pooled.to_csv(a.out + "_pooled.tsv", sep="\t")

    # ---- coverage per nucleus
    cov = pd.concat([pd.read_csv(f, sep="\t") for f in
                     glob.glob(os.path.join(a.dir, "*", "coverage_per_nucleus.tsv"))], ignore_index=True)
    med = cov.groupby("cell_type")[["pos_ge1", "pos_ge2"]].median().round(0)
    med.columns = ["positions_ge1_molecule_median", "positions_ge2_molecules_median"]
    print("\npositions covered per nucleus (median per chromosome, not summed):")
    print(med.to_string())
    med.to_csv(a.out + "_coverage.tsv", sep="\t")

    print(f"\nwrote {a.out}_*.tsv")


if __name__ == "__main__":
    main()
