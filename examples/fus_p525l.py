"""First real-data truth test: FUS P525L knock-in (Smart-seq2, one BAM per cell).

Known truth: P525L cells carry the mutation; control-line cells cannot (used as negative controls).
Usage:
  python examples/fus_p525l.py --meta samples.tsv --star STAR_dir --fasta GRCh38.fa --out out_dir
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd
import pysam

from scmuta.model import SiteFilters, genotype_calls
from scmuta.pileup import ReadFilter, count, discover
from scmuta.pipeline import run_counts

# FUS c.1574C>T p.P525L, GRCh38 (lifted from hg19 chr16:31202752; checked against the reference below)
EXPECTED = ("chr16", 31191431)          # 1-based
WINDOW = ("chr16", 31191300, 31191560)  # 0-based half-open, last exon around codon 525


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", required=True)
    ap.add_argument("--star", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", default="fus_p525l_out")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    # 1. metadata: find the run column (SRR...) and the genotype column (contains P525L)
    meta = pd.read_csv(a.meta, sep="\t")
    print("metadata columns:", list(meta.columns))
    run_col = next(c for c in meta.columns if meta[c].astype(str).str.match(r"^[SED]RR\d+$").mean() > 0.8)
    geno_col = next(c for c in meta.columns if meta[c].astype(str).str.contains("P525L", case=False).any())
    g = meta[geno_col].astype(str)
    group = np.where(g.str.contains("P525L", case=False), "P525L",
                     np.where(g.str.contains(r"KO|knock.?out|-/-", case=False), "KO",
                              np.where(g.str.contains(r"ctrl|control|WT|wild|parental|isogenic", case=False),
                                       "Control", "other")))
    lab = pd.DataFrame({"cell": meta[run_col].astype(str), "genotype": g, "group": group})
    print(lab.groupby(["group", "genotype"]).size().to_string())

    # 2. one BAM per cell; index if needed
    bams = {}
    for c in lab.cell:
        hit = glob.glob(os.path.join(a.star, c, "*.sortedByCoord.out.bam"))
        if hit:
            if not (os.path.exists(hit[0] + ".bai") or os.path.exists(hit[0][:-4] + ".bai")):
                pysam.index(hit[0])
            bams[c] = hit[0]
    lab = lab[lab.cell.isin(bams)]
    print(f"{len(bams)} cells with BAMs:", lab.group.value_counts().to_dict())

    # 3. check the reference base at the expected site
    fa = pysam.FastaFile(a.fasta)
    ref_base = fa.fetch(EXPECTED[0], EXPECTED[1] - 1, EXPECTED[1]).upper()
    ctx = fa.fetch(EXPECTED[0], EXPECTED[1] - 4, EXPECTED[1] + 3).upper()
    print(f"reference at {EXPECTED[0]}:{EXPECTED[1]} = {ref_base} (context {ctx}; expect C, codon CCG>CTG)")

    # 4. scan the window blind to labels: every non-reference allele seen in >= 2 cells
    rf = ReadFilter()
    cand = discover(bams, a.fasta, [WINDOW], rf, None, None, min_alt_reads=3, min_alt_cells=2)
    exp = pd.DataFrame([{"chrom": EXPECTED[0], "pos": EXPECTED[1] - 1, "ref": ref_base,
                         "alt": "T" if ref_base == "C" else "A"}])
    sites = pd.concat([cand[["chrom", "pos", "ref", "alt"]], exp]).drop_duplicates(["chrom", "pos", "alt"])
    sites = sites.reset_index(drop=True)
    A, D, s, cells = count(bams, a.fasta, sites, rf, None, None, dedup="start")

    # per-group alt / depth at each site (what the data says, before any model)
    grp = lab.set_index("cell").group.reindex(cells).to_numpy()
    rows = []
    for k in range(len(s)):
        r = {"pos_1based": int(s.pos[k]) + 1, "ref": s.ref[k], "alt": s.alt[k]}
        for gname in ["P525L", "Control", "KO"]:
            m = grp == gname
            r[f"{gname}_alt/depth"] = f"{A[m, k].sum()}/{D[m, k].sum()}"
            r[f"{gname}_cells_alt/covered"] = f"{(A[m, k] > 0).sum()}/{(D[m, k] > 0).sum()}"
        rows.append(r)
    raw = pd.DataFrame(rows)
    raw.to_csv(os.path.join(a.out, "window_counts.tsv"), sep="\t", index=False)
    print("\nRaw counts in the window (all candidate alleles):")
    print(raw.to_string(index=False))

    # 5. scMuta with control cells as negatives (KO cells kept as cases but scored as non-carriers)
    neg = grp == "Control"
    res = run_counts(A, D, s, cells, negative=neg, filters=SiteFilters(), target_fp=0.01,
                     min_alt_cells=1, out=a.out)
    st = res["sites"]
    print("\nscMuta sites:")
    print(st[["pos", "ref", "alt", "filter", "vaf_pooled", "cells_called", "negative_cells_called"]]
          .assign(pos=st.pos + 1).to_string(index=False))

    # 6. truth check at the P525L site
    k = int(np.flatnonzero((st.pos == EXPECTED[1] - 1).to_numpy())[0])
    post = res["posterior_call"][:, k]
    thr = res["summary"]["threshold"]
    out = pd.DataFrame({"cell": cells, "group": grp, "alt": A[:, k], "depth": D[:, k],
                        "posterior_model": res["posterior"][:, k], "posterior": post})
    # the second (silent) knock-in edit, 4 nt downstream: does each cell agree at both sites?
    k2 = np.flatnonzero(((st.pos == EXPECTED[1] + 3) & (st.alt == "C")).to_numpy())
    if len(k2):
        out["alt_31191435"], out["depth_31191435"] = A[:, k2[0]], D[:, k2[0]]
    out["call"] = genotype_calls(res["posterior"][:, [k]], post[:, None], thr, st["frac_mutant"][k])[:, 0]
    out.to_csv(os.path.join(a.out, "p525l_cells.tsv"), sep="\t", index=False)
    odd = out[((out.group == "Control") & (out.alt > 0)) | ((out.group == "P525L") & (out.call != "mutant"))]
    print("\nCells that disagree with their label or carry any alt read in controls:")
    print(odd.to_string(index=False))
    print(f"\nP525L site: filter={st['filter'][k]}, threshold={thr}")
    print(pd.crosstab(out.group, out.call, margins=True).to_string())
    tab = {gname: out[out.group == gname].call.value_counts().to_dict() for gname in ["P525L", "Control", "KO"]}
    json.dump({"threshold": thr, "site_filter": st["filter"][k], "calls": tab,
               "summary": res["summary"]}, open(os.path.join(a.out, "truth_summary.json"), "w"),
              indent=1, default=float)


if __name__ == "__main__":
    main()
