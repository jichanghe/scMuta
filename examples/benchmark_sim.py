"""Simulation benchmark: scMuta vs a naive per-cell VAF threshold, with and without negative-control cells.

python examples/benchmark_sim.py   (prints a table; 5 random seeds)
"""
import numpy as np
import pandas as pd

from scmuta.pipeline import run_counts
from scmuta.simulate import simulate


def score(called, truth, kinds):
    tp = (called & truth).sum()
    fp = (called & ~truth).sum()
    fn = (~called & truth).sum()
    site_fp = {k: int((called[:, kinds == k].any(0)).sum()) for k in np.unique(kinds) if k != "somatic"}
    return tp / max(tp + fp, 1), tp / max(tp + fn, 1), site_fp


rows = []
for seed in range(5):
    sim = simulate(seed=seed)
    A, D, truth, kinds = sim["A"], sim["D"], sim["truth"], sim["sites"].kind.to_numpy()
    cov_truth = truth & (D >= 2)
    # naive: a cell carries the variant if it has >= 2 alt reads and VAF >= 0.2 (common practice)
    naive = (A >= 2) & (A >= 0.2 * np.maximum(D, 1))
    p, r, sfp = score(naive, cov_truth, kinds)
    rows.append(dict(seed=seed, method="naive VAF>=0.2, >=2 alt", precision=p, recall=r, **sfp))
    for neg, name in ((None, "scMuta, no controls"), (sim["is_ctrl"], "scMuta + negative controls")):
        res = run_counts(A, D, sim["sites"], sim["cells"], negative=neg)
        post = np.nan_to_num(res["posterior_call"], nan=0) >= res["summary"]["threshold"]
        called = post & res["sites"].variant.to_numpy()[None, :]
        p, r, sfp = score(called, cov_truth, kinds)
        rows.append(dict(seed=seed, method=name, precision=p, recall=r, **sfp))
df = pd.DataFrame(rows)
cols = ["precision", "recall", "editing", "germline", "strand", "junction", "error_hotspot", "noise"]
print("Cell-level precision/recall (somatic, cells with >= 2 reads) and number of artefact SITES with any call,"
      " mean of 5 seeds\n")
print(df.groupby("method", sort=False)[cols].mean().round(3).to_string())
