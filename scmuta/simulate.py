"""Count-level simulator with known truth and the main artefacts of single-cell variant calling.

Cells: `n_case` cells that can carry somatic variants (split into clones) and `n_ctrl` negative-control
cells that carry none (another line or donor). Site kinds:
  somatic        - present in one clone of case cells; heterozygous; allelic dropout / monoallelic
                   expression makes some mutant cells show only ref (or only alt) reads
  germline       - heterozygous in every cell (case and control)
  editing        - A>G (T>C) RNA editing at a per-cell level drawn from a beta distribution, every cell
  strand         - alt reads all on one strand (PCR / library artefact), random cells
  junction       - alt reads next to a splice junction (mis-alignment), random cells
  error_hotspot  - site with a raised sequencing error rate, every cell
  noise          - background errors only
Depth per cell-site is negative binomial with many zeros, as in single-cell RNA data.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

KIND_SUBST = {"editing": ("A", "G")}


def simulate(n_case=150, n_ctrl=100, n_clones=3, n_sites=None, depth_mean=6.0, depth_disp=0.5, zero_frac=0.6,
             p_dropout=0.3, base_error=2e-3, seed=0):
    rng = np.random.default_rng(seed)
    n_sites = n_sites or {"somatic": 60, "germline": 40, "editing": 30, "strand": 20, "junction": 20,
                          "error_hotspot": 20, "noise": 80}
    n = n_case + n_ctrl
    is_ctrl = np.r_[np.zeros(n_case, bool), np.ones(n_ctrl, bool)]
    clone = np.r_[rng.integers(0, n_clones, n_case), np.full(n_ctrl, -1)]
    cells = [f"case{i}" for i in range(n_case)] + [f"ctrl{i}" for i in range(n_ctrl)]
    kinds = [k for k, m in n_sites.items() for _ in range(m)]
    S = len(kinds)
    # depth: NB with extra zeros
    p = depth_disp / (depth_disp + depth_mean)
    D = rng.negative_binomial(depth_disp, p, size=(n, S))
    D[rng.random((n, S)) < zero_frac] = 0
    A = np.zeros_like(D)
    truth = np.zeros((n, S), bool)
    feats = np.zeros((S, 8))
    rows = []
    bases = "ACGT"
    for k, kind in enumerate(kinds):
        if kind == "editing":
            ref, alt = "A", "G"
        else:
            ref = bases[rng.integers(4)]
            alt = rng.choice([b for b in bases if b != ref])
        d = D[:, k]
        if kind == "somatic":
            cl = rng.integers(0, n_clones)
            mut = clone == cl
            truth[:, k] = mut
            vaf = np.full(n, 0.5)
            drop = rng.random(n) < p_dropout
            vaf[drop] = rng.choice([0.0, 1.0], size=drop.sum())
            a = np.where(mut, rng.binomial(d, vaf), rng.binomial(d, base_error))
        elif kind == "germline":
            vaf = np.where(rng.random(n) < p_dropout, rng.choice([0.0, 1.0], size=n), 0.5)
            a = rng.binomial(d, vaf)
        elif kind == "editing":
            lvl = rng.beta(0.6, 3.0, size=n)
            a = rng.binomial(d, lvl)
        elif kind in ("strand", "junction"):
            hit = rng.random(n) < 0.25
            a = np.where(hit, rng.binomial(d, 0.35), rng.binomial(d, base_error))
        elif kind == "error_hotspot":
            a = rng.binomial(d, 0.03)
        else:
            a = rng.binomial(d, base_error)
        A[:, k] = a
        alt_tot, ref_tot = a.sum(), (d - a).sum()
        # read features of alt and ref reads
        if kind == "strand":
            af, ar = alt_tot, 0
        else:
            af = rng.binomial(alt_tot, 0.5)
            ar = alt_tot - af
        rf_ = rng.binomial(ref_tot, 0.5)
        near_j = rng.binomial(alt_tot, 0.85 if kind == "junction" else 0.08)
        feats[k] = [af, ar, rf_, ref_tot - rf_, rng.binomial(alt_tot, 0.12), near_j, rng.binomial(alt_tot, 0.03),
                    alt_tot * 36]
        rows.append(("chrS", 1000 * (k + 1), ref, alt, kind))
    sites = pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt", "kind"])
    f = pd.DataFrame(feats, columns=["alt_fwd", "alt_rev", "ref_fwd", "ref_rev", "alt_near_end",
                                     "alt_near_junction", "alt_near_clip", "alt_baseq_sum"])
    sites = pd.concat([sites, f], axis=1)
    sites["alt_reads"] = A.sum(0)
    sites["alt_cells"] = (A > 0).sum(0)
    return {"A": A, "D": D, "sites": sites, "cells": cells, "truth": truth, "is_ctrl": is_ctrl, "clone": clone}
