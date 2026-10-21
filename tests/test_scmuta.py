import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from make_bams import build  # noqa: E402

from scmuta.model import subst_type
from scmuta.pipeline import run_counts
from scmuta.simulate import simulate


def test_subst_type_collapses_strands():
    assert subst_type("A", "G") == "T>C" and subst_type("T", "C") == "T>C" and subst_type("C", "A") == "C>A"


def test_simulation_somatic_recovered_artefacts_rejected():
    sim = simulate(seed=1)
    r = run_counts(sim["A"], sim["D"], sim["sites"], sim["cells"], negative=sim["is_ctrl"])
    v = r["sites"].groupby("kind").variant.sum()
    assert v["somatic"] >= 55
    assert v.drop("somatic").sum() == 0
    post = np.nan_to_num(r["posterior"], nan=0) >= r["summary"]["threshold"]
    called = post & r["sites"].variant.to_numpy()[None, :]
    assert (called & ~sim["truth"]).sum() <= 0.01 * called.sum()     # cell-level precision >= 99%


def test_negative_controls_catch_editing():
    sim = simulate(seed=2)
    with_neg = run_counts(sim["A"], sim["D"], sim["sites"], sim["cells"], negative=sim["is_ctrl"])["sites"]
    no_neg = run_counts(sim["A"], sim["D"], sim["sites"], sim["cells"])["sites"]
    assert with_neg[with_neg.kind == "editing"].variant.sum() == 0
    assert no_neg[no_neg.kind == "editing"].variant.sum() > 10        # without controls editing passes


def test_bam_end_to_end(tmp_path):
    from scmuta.cli import main
    build(str(tmp_path))
    out = tmp_path / "out"
    main(["call", "--plate", str(tmp_path / "cells.tsv"), "--fasta", str(tmp_path / "ref.fa"),
          "--regions", str(tmp_path / "regions.bed"), "--negative", str(tmp_path / "negative.txt"),
          "--out", str(out)])
    s = pd.read_csv(out / "sites.tsv", sep="\t")
    calls = pd.read_csv(out / "cell_calls.tsv", sep="\t")
    v = s[s.variant]
    assert list(zip(v.pos, v.ref, v.alt)) == [(500, "A", "T")]
    assert set(calls.cell) == {f"mut{i}" for i in range(8)}
    assert "splice_junction" in s.loc[s.pos == 999, "filter"].iloc[0]
    assert s.loc[s.pos == 1500, "filter"].iloc[0] == "in_negative_controls"
    assert json.load(open(out / "summary.json"))["variants"] == 1


def test_dedup_counts_molecules(tmp_path):
    from scmuta.pileup import count
    fa, cells = build(str(tmp_path))
    sites = pd.DataFrame({"chrom": ["chrT"], "pos": [500], "ref": ["A"], "alt": ["T"]})
    A1, D1, s1, _ = count(cells, fa, sites)
    A2, D2, s2, _ = count(cells, fa, sites, dedup="start")
    assert D2.sum() <= D1.sum() and A2.sum() <= A1.sum()
    assert {"context", "repeat_extension"} <= set(s1.columns)


def test_not_genotypable_is_not_wild_type():
    from scmuta.model import call_cells, subst_type
    A = np.array([[0], [3], [0]])
    D = np.array([[0], [6], [1]])
    sites = pd.DataFrame({"subst": [subst_type("C", "A")]})
    post, _ = call_cells(sites, A, D, min_depth=2)
    assert np.isnan(post[0, 0]) and np.isnan(post[2, 0]) and post[1, 0] > 0.5


def test_vectorised_em_matches_per_site_loop():
    from scmuta.model import _call_cells_loop, call_cells, site_flags
    sim = simulate(seed=3)
    s = site_flags(sim["sites"], sim["A"], sim["D"])
    p1, s1 = _call_cells_loop(s, sim["A"], sim["D"])
    p2, s2 = call_cells(s, sim["A"], sim["D"])
    assert (np.isnan(p1) == np.isnan(p2)).all()
    f = np.isfinite(p1)
    assert np.abs(p1[f] - p2[f]).max() < 1e-9
    assert np.allclose(s1.frac_mutant, s2.frac_mutant) and np.allclose(s1.mutant_vaf, s2.mutant_vaf)
