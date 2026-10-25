#!/usr/bin/env python
"""Collect the somatic-supported candidate sites from the per-chromosome snrna_pilot runs and describe them.

A candidate passed every site filter, was called in >= 2 nuclei of one donor above the calibrated threshold, and had
>= 5 covered wild-type nuclei of the same donor. This prints what is needed to judge whether they are real:
per-chromosome counts, substitution spectrum (A>G / T>C on the transcribed strand = RNA editing), cell types,
recurrence across donors (a site in several donors is not somatic), and pooled / mutant VAF.

  python examples/snrna_candidates.py --dir ~/snrna_genome --out snrna_candidates.tsv
"""
import argparse
import glob
import os

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", default="snrna_candidates.tsv")
    a = ap.parse_args()

    calls = []
    for f in sorted(glob.glob(os.path.join(os.path.expanduser(a.dir), "chr*", "somatic_calls.tsv"))):
        c = pd.read_csv(f, sep="\t", low_memory=False)
        if "somatic_supported" not in c.columns:
            continue
        # saved as True/False, or as 1.0/0.0 when a donor with no calls made pandas cast the column to float
        c = c[c.somatic_supported.astype(str).str.strip().str.lower().isin(["true", "1", "1.0"])]
        if len(c):
            calls.append(c)
    if not calls:
        print("no somatic-supported calls found")
        return
    x = pd.concat(calls, ignore_index=True)
    print("call columns:", x.columns.tolist())

    key = ["donor", "chrom", "pos"]
    g = x.groupby(key).agg(nuclei=("cell", "nunique"),
                           cell_types=("cell_type", lambda s: ",".join(f"{t}:{n}" for t, n in s.value_counts().items())))
    for col in ("ref", "alt"):
        if col in x.columns:
            g[col] = x.groupby(key)[col].first()
    for col, name in (("alt_reads", "alt_reads_sum"), ("n_alt", "alt_reads_sum"), ("depth", "depth_sum"), ("n_depth", "depth_sum")):
        if col in x.columns and name not in g.columns:
            g[name] = x.groupby(key)[col].sum()
    g = g.reset_index()

    # join site-level stats from each donor's scMuta sites table
    site_cols = ["chrom", "pos", "ref", "alt", "filter", "vaf_pooled", "mutant_vaf", "error_rate", "cells_called",
                 "negative_cells_called", "n_negative_calls"]
    extra = []
    for (d, c), _ in g.groupby(["donor", "chrom"]):
        p = os.path.join(os.path.expanduser(a.dir), c, f"scmuta_{d}", "sites.tsv")
        if not os.path.exists(p):
            continue
        s = pd.read_csv(p, sep="\t", low_memory=False)
        s = s[[k for k in site_cols if k in s.columns]].copy()
        s["donor"] = d
        extra.append(s)
    if extra:
        s = pd.concat(extra, ignore_index=True)
        on = [k for k in ("donor", "chrom", "pos", "ref", "alt") if k in g.columns and k in s.columns]
        s["donor"] = s["donor"].astype(g["donor"].dtype)
        g = g.merge(s.drop_duplicates(on), on=on, how="left")

    if "ref" in g.columns and "alt" in g.columns:
        g["subst"] = g["ref"].astype(str) + ">" + g["alt"].astype(str)

    # same position in more than one donor = inherited or artefact, not somatic
    nd = g.groupby(["chrom", "pos"]).donor.nunique().rename("donors_with_site")
    g = g.merge(nd, left_on=["chrom", "pos"], right_index=True)

    print(f"\n{len(g)} somatic-supported candidate sites")
    print("\nper donor:")
    print(g.groupby("donor").size().to_string())
    print("\nper chromosome:")
    print(g.groupby("chrom").size().sort_index(key=lambda s: s.str.replace("chr", "").str.zfill(2)).to_string())
    if "subst" in g.columns:
        print("\nsubstitutions (A>G and T>C together = the RNA-editing signature):")
        print(g.subst.value_counts().to_string())
        ag = g.subst.isin(["A>G", "T>C"]).mean()
        print(f"A>G + T>C share: {ag:.0%}  (random expectation ~17%; >50% points to editing)")
    print("\nsites present in more than one donor:", int((g.donors_with_site > 1).sum()))
    print("\nprimary cell type of each site:")
    print(g.cell_types.str.split(":").str[0].value_counts().to_string())
    print("\nnuclei per site:")
    print(g.nuclei.describe().round(1).to_string())
    show = [c for c in ["donor", "chrom", "pos", "subst", "nuclei", "cell_types", "vaf_pooled", "mutant_vaf",
                        "donors_with_site", "filter"] if c in g.columns]
    print("\ntop 40 by nuclei:")
    out = g.sort_values("nuclei", ascending=False)
    print(out[show].assign(pos=out.pos + 1).head(40).to_string(index=False))
    out.to_csv(a.out, sep="\t", index=False)
    print(f"\nwrote {a.out} (pos is 0-based, as in scMuta tables)")


if __name__ == "__main__":
    main()
