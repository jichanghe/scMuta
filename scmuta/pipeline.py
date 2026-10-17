"""End-to-end run on count matrices: flag sites, call cells, calibrate on negative controls, check reproducibility."""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

from .model import SiteFilters, calibrate_threshold, call_cells, site_flags, split_check
from .util import log


def run_counts(A, D, sites, cells, negative=None, filters: SiteFilters | None = None, target_fp=0.001,
               threshold=None, min_alt_cells=2, neg_p=1e-3, seed=0, out=None):
    """A, D: cells x sites alt and depth counts. negative: boolean mask of negative-control cells (or None).

    Returns dict with site table, long call table, calibration table and summary."""
    s = site_flags(sites, A, D, filters)
    post, s = call_cells(s, A, D)
    usable = (s["filter"] == "PASS").to_numpy()
    calib = None
    case = np.ones(len(cells), bool)
    if negative is not None and negative.any():
        negative = np.asarray(negative, bool)
        case = ~negative
        # cross-fitting on the negative controls: split them into two folds; each fold's false-call rate is
        # measured at sites screened with the OTHER fold only, so the reported rate is not estimated on the
        # cells that removed the artefact sites. The final site screen uses all negative controls (most power).
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
        tabs = []
        for k in (0, 1):
            other_bad = screen(folds[1 - k])[0]
            held = np.zeros(len(cells), bool)
            held[folds[k]] = True
            tabs.append(calibrate_threshold(post, held, usable & ~other_bad, target_fp)[1])
        calib = tabs[0].copy()
        calib["negative_calls"] = tabs[0].negative_calls + tabs[1].negative_calls
        calib["negative_pairs"] = tabs[0].negative_pairs + tabs[1].negative_pairs
        calib["negative_fp_rate"] = calib.negative_calls / calib.negative_pairs.clip(lower=1)
        ok = calib[calib.negative_fp_rate <= target_fp]
        thr = float(ok.threshold.min()) if len(ok) else float(calib.threshold.max())
        s.loc[bad_all & usable, "filter"] = "in_negative_controls"
        usable = usable & ~bad_all
    else:
        thr = 0.9
    if threshold is not None:
        thr = threshold
    called = (np.nan_to_num(post, nan=0) >= thr)
    s["cells_called"] = called[case].sum(0)
    if negative is not None and negative.any():
        s["negative_cells_called"] = called[negative].sum(0)
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
        row = calib[calib.threshold == thr]
        summary["held_out_negative_fp_rate"] = float(row.negative_fp_rate.iloc[0]) if len(row) else None
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
    return {"sites": s, "calls": calls, "posterior": post, "calibration": calib, "summary": summary}
