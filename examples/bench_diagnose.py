"""Why does scMuta miss carrier pairs, and what are the control-cell calls? (reads cellline_benchmark.py outputs)

  python examples/bench_diagnose.py --bench bench_cellline
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cellline_benchmark as cb  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", default="bench_cellline")
    a = ap.parse_args()
    t = pd.read_csv(os.path.join(a.bench, "cells_split.tsv"), sep="\t", dtype=str)
    z = np.load(os.path.join(a.bench, "counts.npz"), allow_pickle=True)
    A, D = sparse.csr_matrix(z["A"].item()), sparse.csr_matrix(z["D"].item())
    sites = pd.read_csv(os.path.join(a.bench, "sites_counted.tsv"), sep="\t")
    line = t["cell line"].to_numpy()
    truth = (t.half == "truth").to_numpy()
    lab, _ = cb.truth_labels(A, D, line, truth)
    tsites = np.flatnonzero(sites.found_in.fillna("").str.contains("test").to_numpy())
    st = pd.read_csv(os.path.join(a.bench, "scmuta", "sites.tsv"), sep="\t")
    st["truth"] = lab[tsites]
    tc = np.flatnonzero(~truth & np.isin(line, cb.CASE_LINES + cb.CONTROL_LINES))
    lt = line[tc]
    At = A[tc][:, tsites].tocsc()
    Dt = D[tc][:, tsites].tocsc()

    # 1. carrier pairs (cell of the owning line, >= 2 alt reads: callable) by the site's scMuta filter
    print("\n1. callable carrier pairs (>= 2 alt reads) by scMuta site filter")
    rows = []
    for L in cb.CASE_LINES:
        cols = np.flatnonzero((st.truth == f"private_{L}").to_numpy())
        m = lt == L
        n2 = np.asarray((At[m][:, cols] >= 2).sum(0)).ravel()
        rows.append(pd.DataFrame({"filter": st["filter"].to_numpy()[cols], "variant": st["variant"].to_numpy()[cols],
                                  "pairs": n2}))
    r = pd.concat(rows)
    r["first_filter"] = r["filter"].astype(str).str.split(",").str[0]
    g = r.groupby(["variant", "first_filter"]).pairs.sum().sort_values(ascending=False)
    print((g / g.sum()).round(3).to_string())
    print("all filters (a site can have several):")
    ex = r.assign(f=r["filter"].astype(str).str.split(",")).explode("f")
    print(ex.groupby("f").pairs.sum().sort_values(ascending=False).head(12).to_string())
    for col in ("vaf_pooled", "frac_mutant", "mutant_vaf", "cells_called"):
        if col in st.columns:
            sub = st[st.truth.isin([f"private_{L}" for L in cb.CASE_LINES])]
            print(f"private sites {col}: by PASS / not:",
                  sub.groupby(sub["filter"] == "PASS")[col].median().round(3).to_dict())

    # 2. control cells at sites absent from every line: per-cell noise, VAF, recurrence
    print("\n2. control (HSPC) test cells at 'absent' sites, raw counts")
    absent = np.flatnonzero((st.truth == "absent").to_numpy() & (st["filter"] == "PASS").to_numpy())
    for L in cb.CASE_LINES + cb.CONTROL_LINES:
        m = np.flatnonzero(lt == L)
        a_ = At[m][:, absent].toarray()
        d_ = Dt[m][:, absent].toarray()
        hit = (a_ >= 2) & (a_ >= 0.2 * d_)
        per_cell = hit.sum(1) / np.maximum((d_ > 0).sum(1), 1)
        vaf = (a_[hit] / d_[hit]) if hit.any() else np.array([np.nan])
        rec = hit.sum(0)
        print(f"{L}: cells {len(m)}; '>=2 alt & VAF>=0.2' per covered pair: median {np.median(per_cell):.2e}, "
              f"max {per_cell.max():.2e}; hit VAF median {np.nanmedian(vaf):.2f}; "
              f"sites hit in 1 cell {int((rec == 1).sum())}, in >=2 cells {int((rec >= 2).sum())}")
        if L in cb.CONTROL_LINES:
            print("   per-cell rate, sorted:", np.round(np.sort(per_cell)[::-1] * 1e4, 2).tolist(), "(x1e-4)")
            subst = st.iloc[absent].apply(lambda x: f"{x.ref}>{x.alt}", axis=1).to_numpy()
            s_hit = pd.Series(np.repeat(subst[None, :], len(m), 0)[hit]).value_counts()
            print("   substitutions of hits:", s_hit.head(8).to_dict())


if __name__ == "__main__":
    main()
