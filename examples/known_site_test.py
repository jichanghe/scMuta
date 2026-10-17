"""Truth test at known mutation sites in plate data (one BAM per cell) with labelled cell groups.

Carrier groups carry the mutation; negative groups cannot (they are scMuta's negative controls).
Around each known site, every non-reference allele is also scanned blind to labels.

Example (TARGET-seq cell lines, GSE105451):
  python examples/known_site_test.py --runs metadata/runs.tsv --star STAR \
      --fasta GRCh38.fa --site JAK2_V617F:chr9:5073770:G:T --carriers SET2 --negatives JURKAT,HSPC --out out/
"""
import argparse
import glob
import json
import os
import re

import numpy as np
import pandas as pd
import pysam

from scmuta.model import SiteFilters
from scmuta.pileup import ReadFilter, count, discover
from scmuta.pipeline import run_counts


def groups_from_runs(path):
    """SRA run table (pysradb --detailed): run_accession + group parsed from the sample title,
    e.g. 'GSM2826274: JURKAT_CLPL2_10B; Homo sapiens; RNA-Seq' -> JURKAT."""
    m = pd.read_csv(path, sep="\t", dtype=str)
    title_col = next(c for c in ["experiment_title", "sample_title", "experiment_alias"] if c in m.columns)
    def grp(t):
        t = str(t).split(": ", 1)[-1].split(";")[0]
        return re.split(r"[_\s-]", t)[0].upper()
    return pd.DataFrame({"cell": m["run_accession"], "group": m[title_col].map(grp),
                         "title": m[title_col]}).drop_duplicates("cell")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True, help="SRA run table (pysradb metadata --detailed)")
    ap.add_argument("--star", required=True, help="folder with <run>/Aligned.sortedByCoord.out.bam")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--site", action="append", required=True, help="NAME:chrom:pos(1-based):ref:alt, repeatable")
    ap.add_argument("--carriers", required=True, help="comma-separated carrier groups")
    ap.add_argument("--negatives", required=True, help="comma-separated negative-control groups")
    ap.add_argument("--window", type=int, default=150)
    ap.add_argument("--out", default="known_site_out")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    carriers = [x.strip().upper() for x in a.carriers.split(",")]
    negatives = [x.strip().upper() for x in a.negatives.split(",")]

    lab = groups_from_runs(a.runs)
    bams = {}
    for c in lab.cell:
        hit = glob.glob(os.path.join(a.star, c, "*.sortedByCoord.out.bam"))
        if hit and os.path.exists(hit[0] + ".bai"):
            bams[c] = hit[0]
    lab = lab[lab.cell.isin(bams)]
    print(f"{len(bams)} cells with indexed BAMs:", lab.group.value_counts().to_dict())
    if not set(carriers) & set(lab.group) or not set(negatives) & set(lab.group):
        print("warning: carriers or negatives missing among aligned cells so far")

    fa = pysam.FastaFile(a.fasta)
    known, windows = [], []
    for spec in a.site:
        name, chrom, pos, ref, alt = spec.split(":")
        pos = int(pos)
        rb = fa.fetch(chrom, pos - 1, pos).upper()
        ctx = fa.fetch(chrom, pos - 6, pos + 5).upper()
        print(f"{name}: reference at {chrom}:{pos} = {rb} (expected {ref}); context {ctx[:5]}[{ctx[5]}]{ctx[6:]}")
        if rb != ref.upper():
            print(f"  WARNING: reference base does not match; check the coordinate for {name}")
        known.append({"name": name, "chrom": chrom, "pos": pos - 1, "ref": rb, "alt": alt.upper()})
        windows.append((chrom, pos - 1 - a.window, pos + a.window))
    known = pd.DataFrame(known)

    rf = ReadFilter()
    cand = discover(bams, a.fasta, windows, rf, None, None, min_alt_reads=3, min_alt_cells=2)
    sites = pd.concat([cand[["chrom", "pos", "ref", "alt"]], known[["chrom", "pos", "ref", "alt"]]])
    sites = sites.drop_duplicates(["chrom", "pos", "alt"]).sort_values(["chrom", "pos"]).reset_index(drop=True)
    A, D, s, cells = count(bams, a.fasta, sites, rf, None, None, dedup="start")
    grp = lab.set_index("cell").group.reindex(cells).to_numpy()
    groups = carriers + negatives + sorted(set(lab.group) - set(carriers) - set(negatives))

    rows = []
    for k in range(len(s)):
        r = {"pos_1based": int(s.pos[k]) + 1, "ref": s.ref[k], "alt": s.alt[k]}
        for g in groups:
            m = grp == g
            r[f"{g} alt/depth"] = f"{A[m, k].sum()}/{D[m, k].sum()}"
            r[f"{g} cells alt/cov"] = f"{(A[m, k] > 0).sum()}/{(D[m, k] > 0).sum()}"
        rows.append(r)
    raw = pd.DataFrame(rows)
    raw.to_csv(os.path.join(a.out, "window_counts.tsv"), sep="\t", index=False)
    print("\nRaw counts (all candidate alleles, found blind to labels, plus the known sites):")
    print(raw.to_string(index=False))

    neg = np.isin(grp, negatives)
    res = run_counts(A, D, s, cells, negative=neg, filters=SiteFilters(), target_fp=0.01,
                     min_alt_cells=1, out=a.out)
    st = res["sites"]
    print("\nscMuta sites:")
    print(st[["pos", "ref", "alt", "filter", "vaf_pooled", "cells_called", "negative_cells_called"]]
          .assign(pos=st.pos + 1).to_string(index=False))
    thr = res["summary"]["threshold"]
    print(f"\nposterior threshold: {thr}")

    summary = {"threshold": thr, "sites": {}}
    for _, kn in known.iterrows():
        k = int(np.flatnonzero(((st.pos == kn.pos) & (st.alt == kn.alt)).to_numpy())[0])
        post = res["posterior_call"][:, k]
        out = pd.DataFrame({"cell": cells, "group": grp, "alt": A[:, k], "depth": D[:, k],
                            "posterior_model": res["posterior"][:, k], "posterior": post})
        out["call"] = np.where(np.isnan(post), "not_genotypable", np.where(post >= thr, "mutant", "wild_type"))
        out.to_csv(os.path.join(a.out, f"{kn['name']}_cells.tsv"), sep="\t", index=False)
        print(f"\n=== {kn['name']} ({kn.chrom}:{kn.pos + 1} {kn.ref}>{kn.alt}), filter={st['filter'][k]}")
        print(pd.crosstab(out.group, out.call, margins=True).to_string())
        odd = out[(np.isin(out.group, negatives) & (out.call == "mutant")) |
                  (np.isin(out.group, carriers) & (out.call == "wild_type"))]
        if len(odd):
            print("cells disagreeing with their label:")
            print(odd.to_string(index=False))
        summary["sites"][kn["name"]] = {"filter": st["filter"][k],
                                        "calls": {g: out[out.group == g].call.value_counts().to_dict() for g in groups}}
    summary["scmuta"] = res["summary"]
    json.dump(summary, open(os.path.join(a.out, "truth_summary.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
