"""End-to-end run on count matrices: flag sites, call cells, calibrate on negative controls, check reproducibility."""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

from .model import (SiteFilters, calibrate_threshold, call_cells, callable_posterior, site_flags, site_low_vaf_counts,
                    split_check, subst_type, substitution_error_from_sums)
from .util import log


GRID = np.array([0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999])


def loo_negative_calls(pc, A, D, neg_idx, usable, bg, neg_p=1e-3, grid=GRID, case=None, min_alt_cells=1):
    """Leave-one-out calibration on negative-control cells.

    For each negative cell, a site counts towards its false-call rate only if the OTHER negative cells do not
    already carry the allele (binomial excess over the background error with >= 2 alt reads). This removes the
    controls' own germline variants and shared artefacts from the error estimate (the final site screen removes
    them from the calls) without using the cell's own reads, so the estimate is not circular.
    A control cell's call counts only at sites that would be reported at that threshold (called in >=
    min_alt_cells case cells), so the estimate is the false-call rate of the reported calls; single-cell noise at
    sites that no case cell shares is never reported and is not counted.
    Returns (calls per threshold, covered pairs)."""
    from scipy import stats
    an = np.asarray(A[neg_idx], float)
    dn = np.asarray(D[neg_idx], float)
    a_oth = an.sum(0)[None, :] - an
    d_oth = dn.sum(0)[None, :] - dn
    with np.errstate(invalid="ignore"):
        p = stats.binom.sf(a_oth - 1, d_oth, np.broadcast_to(bg, a_oth.shape))
    bad = (p < neg_p) & (a_oth >= 2)
    sub = np.array(pc[neg_idx], float)
    sub[:, ~usable] = np.nan
    sub[bad] = np.nan
    cov = np.isfinite(sub)
    sub0 = np.nan_to_num(sub, nan=0)
    if case is None or min_alt_cells <= 1:
        return np.array([(sub0 >= t).sum() for t in grid], float), float(cov.sum())
    pcc = np.nan_to_num(np.asarray(pc[np.asarray(case, bool)], float), nan=0)
    out = []
    for t in grid:
        reported = (pcc >= t).sum(0) >= min_alt_cells
        out.append((sub0[:, reported] >= t).sum())
    return np.array(out, float), float(cov.sum())


def run_counts(A, D, sites, cells, negative=None, filters: SiteFilters | None = None, target_fp=0.001,
               threshold=None, min_alt_cells=2, neg_p=1e-3, seed=0, out=None, min_cell_alt=2, min_threshold=0.9):
    """A, D: cells x sites alt and depth counts. negative: boolean mask of negative-control cells (or None).

    Returns dict with site table, long call table, calibration table and summary."""
    s = site_flags(sites, A, D, filters)
    post, s = call_cells(s, A, D)
    post_call = callable_posterior(post, A, min_cell_alt)
    usable = (s["filter"] == "PASS").to_numpy()
    calib = None
    case = np.ones(len(cells), bool)
    if negative is not None and negative.any():
        negative = np.asarray(negative, bool)
        case = ~negative
        # leave-one-out on the negative controls: each control cell's false-call rate is measured at sites screened
        # with the OTHER control cells only (see loo_negative_calls), so the rate is not estimated on the reads
        # that removed the site. The final site screen uses all negative controls (most power).
        from scipy import stats
        rng = np.random.default_rng(seed)
        neg_idx = np.flatnonzero(negative)
        rng.shuffle(neg_idx)
        folds = [neg_idx[: len(neg_idx) // 2], neg_idx[len(neg_idx) // 2:]]
        bg = 2 * s["background_error"].to_numpy()

        def screen(idx):
            a_n, d_n = A[idx].sum(0), D[idx].sum(0)
            p = stats.binom.sf(a_n - 1, d_n, bg)
            return (p < neg_p) & (a_n >= 2), a_n, d_n, p

        bad_all, a_n, d_n, p_neg = screen(neg_idx)
        s["negative_alt_reads"], s["negative_depth"], s["negative_p"] = a_n, d_n, p_neg
        nc, npairs = loo_negative_calls(post_call, A, D, neg_idx, usable, bg, neg_p, case=case,
                                        min_alt_cells=min_alt_cells)
        calib = pd.DataFrame({"threshold": GRID, "negative_calls": nc.astype(int), "negative_pairs": int(npairs),
                              "negative_fp_rate": nc / max(npairs, 1)})
        ok = calib[calib.negative_fp_rate <= target_fp]
        thr_calibrated = float(ok.threshold.min()) if len(ok) else float(calib.threshold.max())
        # the calibrated threshold can only raise the bar: when negatives give no calls at all, the lowest grid
        # value passes, which is too permissive for cells outside the controls
        thr = max(thr_calibrated, min_threshold)
        s.loc[bad_all & usable, "filter"] = "in_negative_controls"
        usable = usable & ~bad_all
    else:
        thr = 0.9
    if threshold is not None:
        thr = threshold
    called = (np.nan_to_num(post_call, nan=0) >= thr)
    s["cells_called"] = called[case].sum(0)
    if negative is not None and negative.any():
        s["negative_cells_called"] = called[negative].sum(0)
        s["negative_cells_covered"] = (D[negative] >= 2).sum(0)
    # a site needs calls in >= min_alt_cells cells of interest to be reported as a variant
    s["variant"] = usable & (s["cells_called"] >= min_alt_cells)
    rep = split_check(A, D, s, s["variant"].to_numpy(), thr)
    ii, kk = np.nonzero(called & s["variant"].to_numpy()[None, :] & case[:, None])
    calls = pd.DataFrame({"cell": np.asarray(cells)[ii], "chrom": s.chrom.to_numpy()[kk], "pos": s.pos.to_numpy()[kk],
                          "ref": s.ref.to_numpy()[kk], "alt": s.alt.to_numpy()[kk],
                          "alt_reads": A[ii, kk], "depth": D[ii, kk], "posterior": post[ii, kk]})
    summary = {"cells": len(cells), "candidate_sites": int(len(s)), "pass_filters": int(usable.sum()),
               "variants": int(s["variant"].sum()), "cell_calls": int(len(calls)), "threshold": thr,
               "filters": s["filter"].str.split(",").explode().value_counts().to_dict(),
               "split_check": rep}
    if calib is not None:
        summary["negative_controls"] = int(negative.sum())
        summary["threshold_calibrated"] = thr_calibrated
        row = calib[calib.threshold == thr]
        summary["held_out_negative_fp_rate"] = float(row.negative_fp_rate.iloc[0]) if len(row) else None
        # the pooled rate can hide a site where negatives are called; list those sites explicitly
        bad = s[s["variant"] & (s["negative_cells_called"] > 0)]
        summary["variant_sites_with_negative_calls"] = int(len(bad))
    log(f"variants: {summary['variants']} (of {len(s)} candidates; {usable.sum()} pass filters); "
        f"{len(calls)} cell calls at posterior >= {thr}")
    if out:
        os.makedirs(out, exist_ok=True)
        s.to_csv(os.path.join(out, "sites.tsv"), sep="\t", index=False)
        calls.to_csv(os.path.join(out, "cell_calls.tsv"), sep="\t", index=False)
        if calib is not None:
            calib.to_csv(os.path.join(out, "calibration.tsv"), sep="\t", index=False)
        with open(os.path.join(out, "summary.json"), "w") as fh:
            json.dump(summary, fh, indent=2, default=float)
    return {"sites": s, "calls": calls, "posterior": post, "posterior_call": post_call, "calibration": calib, "summary": summary}


def run_counts_chunked(A, D, sites, cells, negative, filters: SiteFilters | None = None, target_fp=0.001,
                       min_alt_cells=2, neg_p=1e-3, seed=0, out=None, min_cell_alt=2, min_threshold=0.9,
                       chunk=1500, store_min=None, max_block_bytes=2e8):
    """Memory-light version of run_counts for large sparse data (e.g. 10x snRNA: 10^4-10^5 cells x 10^4-10^5
    sites). A, D: scipy sparse (cells x sites). Sites are processed in blocks; the negative-control calibration is
    pooled over all blocks before the threshold is chosen, so it is the same procedure as run_counts.

    Memory is bounded two ways, because the posteriors of every covered cell-site pair over a whole chromosome
    do not fit for large data:
      store_min       posteriors below this are not kept (they can never be called: the chosen threshold is at
                      least min_threshold). Default min_threshold. Pass 0.0 to keep them all, which callers
                      that need the full posterior matrix - a threshold sweep, wild-type calls - must do.
      max_block_bytes caps the dense working array per block, reducing `chunk` when there are many cells."""
    from scipy import sparse, stats
    A = sparse.csc_matrix(A)
    D = sparse.csc_matrix(D)
    negative = np.asarray(negative, bool)
    case = ~negative
    rng = np.random.default_rng(seed)
    neg_idx = np.flatnonzero(negative)
    rng.shuffle(neg_idx)
    folds = [neg_idx[: len(neg_idx) // 2], neg_idx[len(neg_idx) // 2:]]
    held = [np.isin(np.arange(len(cells)), f) for f in folds]
    grid = np.array([0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999])
    neg_calls = np.zeros((2, len(grid)))
    neg_pairs = np.zeros(2)
    blocks = []
    sites = sites.reset_index(drop=True)
    store_min = min_threshold if store_min is None else store_min
    chunk = max(100, min(chunk, int(max_block_bytes / max(8 * A.shape[0], 1))))
    # background error per substitution class from ALL sites (so the block size cannot change results)
    a_low, d_low = np.zeros(len(sites)), np.zeros(len(sites))
    for b0 in range(0, len(sites), chunk):
        sl = slice(b0, min(len(sites), b0 + chunk))
        a_low[sl], d_low[sl] = site_low_vaf_counts(A[:, sl].toarray(), D[:, sl].toarray())
    subst = [subst_type(r, a_) for r, a_ in zip(sites.ref, sites.alt)]
    sub_err = substitution_error_from_sums(subst, a_low, d_low)
    for b0 in range(0, len(sites), chunk):
        sl = slice(b0, min(len(sites), b0 + chunk))
        a = A[:, sl].toarray()
        d = D[:, sl].toarray()
        s = site_flags(sites.iloc[sl].reset_index(drop=True), a, d, filters)
        post, s = call_cells(s, a, d, sub_err=sub_err)
        pc = callable_posterior(post, a, min_cell_alt)
        usable = (s["filter"] == "PASS").to_numpy()
        bg = 2 * s["background_error"].to_numpy()

        def screen(idx):
            a_n, d_n = a[idx].sum(0), d[idx].sum(0)
            p = stats.binom.sf(a_n - 1, d_n, bg)
            return (p < neg_p) & (a_n >= 2), a_n, d_n, p

        bad_all, a_n, d_n, p_neg = screen(neg_idx)
        s["negative_alt_reads"], s["negative_depth"], s["negative_p"] = a_n, d_n, p_neg
        nc, npairs = loo_negative_calls(pc, a, d, neg_idx, usable, bg, neg_p, grid, case=case,
                                        min_alt_cells=min_alt_cells)
        neg_calls[0] += nc
        neg_pairs[0] += npairs
        s.loc[bad_all & usable, "filter"] = "in_negative_controls"
        # keep posteriors sparsely: covered pairs that could still be called at the chosen threshold
        keep = np.isfinite(pc) & (np.nan_to_num(pc, nan=0) >= store_min)
        rows, cols = np.nonzero(keep)
        blocks.append((sl, s, sparse.csc_matrix((pc[rows, cols], (rows, cols)), shape=pc.shape),
                       sparse.csc_matrix((post[rows, cols] + 1e-12, (rows, cols)), shape=post.shape)))
        del a, d, post, pc, keep, rows, cols
        log(f"sites {sl.start}-{sl.stop} of {len(sites)}")
    rate = neg_calls.sum(0) / max(neg_pairs.sum(), 1)
    calib = pd.DataFrame({"threshold": grid, "negative_calls": neg_calls.sum(0).astype(int),
                          "negative_pairs": int(neg_pairs.sum()), "negative_fp_rate": rate})
    ok = calib[calib.negative_fp_rate <= target_fp]
    thr_cal = float(ok.threshold.min()) if len(ok) else float(grid.max())
    thr = max(thr_cal, min_threshold)
    site_tabs, call_tabs = [], []
    for sl, s, pcs, _ in blocks:
        called = (pcs >= thr).toarray()
        s["cells_called"] = called[case].sum(0)
        s["negative_cells_called"] = called[negative].sum(0)
        s["variant"] = (s["filter"] == "PASS").to_numpy() & (s["cells_called"] >= min_alt_cells)
        ii, kk = np.nonzero(called & s["variant"].to_numpy()[None, :] & case[:, None])
        if len(ii):
            a_ = A[:, sl]
            d_ = D[:, sl]
            call_tabs.append(pd.DataFrame({
                "cell": np.asarray(cells)[ii], "chrom": s.chrom.to_numpy()[kk], "pos": s.pos.to_numpy()[kk],
                "ref": s.ref.to_numpy()[kk], "alt": s.alt.to_numpy()[kk],
                "alt_reads": np.asarray(a_[ii, kk]).ravel(), "depth": np.asarray(d_[ii, kk]).ravel(),
                "posterior": np.asarray(pcs[ii, kk]).ravel()}))
        site_tabs.append(s)
    st = pd.concat(site_tabs, ignore_index=True)
    calls = pd.concat(call_tabs, ignore_index=True) if call_tabs else pd.DataFrame(
        columns=["cell", "chrom", "pos", "ref", "alt", "alt_reads", "depth", "posterior"])
    summary = {"cells": len(cells), "candidate_sites": int(len(st)),
               "pass_filters": int((st["filter"] == "PASS").sum()), "variants": int(st["variant"].sum()),
               "cell_calls": int(len(calls)), "threshold": thr, "threshold_calibrated": thr_cal,
               "negative_controls": int(negative.sum()),
               "held_out_negative_fp_rate": float(rate[grid == thr][0]) if (grid == thr).any() else None,
               "variant_sites_with_negative_calls": int((st["variant"] & (st["negative_cells_called"] > 0)).sum()),
               "filters": st["filter"].str.split(",").explode().value_counts().to_dict()}
    log(f"variants: {summary['variants']} (of {len(st)} candidates); {len(calls)} cell calls at posterior >= {thr}")
    if out:
        os.makedirs(out, exist_ok=True)
        st.to_csv(os.path.join(out, "sites.tsv"), sep="\t", index=False)
        calls.to_csv(os.path.join(out, "cell_calls.tsv"), sep="\t", index=False)
        calib.to_csv(os.path.join(out, "calibration.tsv"), sep="\t", index=False)
        with open(os.path.join(out, "summary.json"), "w") as fh:
            json.dump(summary, fh, indent=2, default=float)
    # posteriors of covered cell-site pairs (sparse; 0 = not covered), for genotype states downstream
    post_call_sp = sparse.hstack([b[2] for b in blocks]).tocsc() if blocks else None
    post_sp = sparse.hstack([b[3] for b in blocks]).tocsc() if blocks else None
    return {"sites": st, "calls": calls, "calibration": calib, "summary": summary,
            "posterior_call_sparse": post_call_sp, "posterior_sparse": post_sp}
