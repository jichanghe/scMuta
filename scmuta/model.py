"""Site-level artefact flags, per-cell genotype posteriors, and calibration against negative-control cells.

Model per candidate site s and cell i (a = alt reads, d = depth):
    wild type : a ~ Binomial(d, eps_s)            eps_s = site error, shrunk toward the error of its
                                                  substitution type across all sites
    mutant    : a ~ BetaBinomial(d, alpha, beta)  mean VAF ~ 0.5 by default, wide enough for allelic
                                                  imbalance and partial dropout in single cells
    P(mutant) = pi_s, fitted by EM per site
The per-cell posterior P(mutant | a, d) is the output. Cells with d = 0 get no call.

Calibration: when some cells are known not to carry the variants of interest (negative controls: another
cell line, another donor, cells before an engineering step), the posterior threshold is chosen so the
negative-control false-call rate stays below a target, and the report shows how many calls negatives get.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import betaln, gammaln

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def subst_type(ref: str, alt: str) -> str:
    """Strand-collapsed substitution class, e.g. A>G and T>C -> T>C (pyrimidine reference)."""
    if ref in "AG":
        ref, alt = ref.translate(COMPLEMENT), alt.translate(COMPLEMENT)
    return f"{ref}>{alt}"


@dataclass
class SiteFilters:
    stranded_protocol: bool = False    # 10x 3'/5': reads come from one strand; strand bias is meaningless
    max_strand_bias_p: float = 1e-3    # Fisher p for alt vs ref strand distribution
    min_strand_shift: float = 0.25     # ... and the forward fraction of alt reads must differ from ref by this
    max_frac_near_end: float = 0.6     # alt reads within 10 nt of a read end
    max_frac_near_junction: float = 0.5  # alt reads within 10 nt of a splice junction in the read
    max_frac_near_clip: float = 0.5
    min_mean_baseq: float = 25.0
    germline_vaf: float = 0.3          # pooled VAF above this and present in most covered cells -> germline
    germline_cell_frac: float = 0.6
    editing_bed: str | None = None     # known RNA-editing sites (REDIportal etc.); A>G/T>C flagged if listed
    flag_all_ag: bool = False          # flag every A>G / T>C as possible editing when no list is given
    editing_cluster_window: int = 0    # >0: flag A>G/T>C candidates with another within this many nt


def site_flags(sites: pd.DataFrame, A: np.ndarray, D: np.ndarray, f: SiteFilters | None = None) -> pd.DataFrame:
    f = f or SiteFilters()
    s = sites.copy()
    alt_tot = A.sum(0)
    dep_tot = D.sum(0)
    cov = (D > 0).sum(0)
    with_alt = (A > 0).sum(0)
    s["depth"] = dep_tot
    s["alt_total"] = alt_tot
    s["vaf_pooled"] = np.where(dep_tot > 0, alt_tot / np.maximum(dep_tot, 1), np.nan)
    s["cells_covered"] = cov
    s["cells_alt"] = with_alt
    s["subst"] = [subst_type(r, a) for r, a in zip(s.ref, s.alt)]
    n_alt = np.maximum(alt_tot, 1)
    flags = []
    clustered = editing_clusters(s, f.editing_cluster_window) if f.editing_cluster_window else np.zeros(len(s), bool)
    edit_sites = set()
    if f.editing_bed:
        bed = pd.read_csv(f.editing_bed, sep="\t", header=None, usecols=[0, 1], comment="#")
        edit_sites = set(zip(bed[0].astype(str), bed[1].astype(int)))
    for k in range(len(s)):
        fl = []
        if not f.stranded_protocol and alt_tot[k] >= 4:
            p = stats.fisher_exact([[s.alt_fwd[k], s.alt_rev[k]], [s.ref_fwd[k] + 1, s.ref_rev[k] + 1]])[1]
            # a p-value alone flags real variants at high depth for small imbalances; also require a large shift
            f_alt = s.alt_fwd[k] / max(s.alt_fwd[k] + s.alt_rev[k], 1)
            f_ref = (s.ref_fwd[k] + 1) / (s.ref_fwd[k] + s.ref_rev[k] + 2)
            if s.ref_fwd[k] + s.ref_rev[k] < 10:   # homozygous / few ref reads: compare with an even split
                f_ref = 0.5
            one_sided = min(s.alt_fwd[k], s.alt_rev[k]) == 0 and alt_tot[k] >= 6
            if (p < f.max_strand_bias_p and abs(f_alt - f_ref) > f.min_strand_shift) or one_sided:
                fl.append("strand_bias")
        if s.alt_near_end[k] / n_alt[k] > f.max_frac_near_end:
            fl.append("read_end")
        if s.alt_near_junction[k] / n_alt[k] > f.max_frac_near_junction:
            fl.append("splice_junction")
        if s.alt_near_clip[k] / n_alt[k] > f.max_frac_near_clip:
            fl.append("soft_clip")
        if alt_tot[k] and s.alt_baseq_sum[k] / alt_tot[k] < f.min_mean_baseq:
            fl.append("low_baseq")
        if "repeat_extension" in s.columns and s.repeat_extension[k]:
            fl.append("repeat_extension")
        if s.vaf_pooled[k] > f.germline_vaf and cov[k] > 0 and with_alt[k] / cov[k] > f.germline_cell_frac:
            fl.append("germline_like")
        if s.subst[k] == "T>C":
            if (s.chrom[k], int(s.pos[k]) + 1) in edit_sites or (s.chrom[k], int(s.pos[k])) in edit_sites:
                fl.append("rna_editing_known")
            elif f.flag_all_ag:
                fl.append("possible_editing")
        if clustered[k]:
            fl.append("editing_cluster")
        flags.append(",".join(fl) if fl else "PASS")
    s["filter"] = flags
    return s


def editing_clusters(sites: pd.DataFrame, window: int = 50) -> np.ndarray:
    """A>G / T>C candidates with another A>G / T>C candidate within `window` nt: the clustered pattern of A-to-I
    RNA editing (hyper-edited Alu and other dsRNA). Somatic SNVs rarely come in tight same-type clusters."""
    st = np.array([subst_type(r, a) for r, a in zip(sites.ref, sites.alt)])
    ag = st == "T>C"
    out = np.zeros(len(sites), bool)
    for c in pd.unique(sites.chrom):
        idx = np.flatnonzero((sites.chrom.to_numpy() == c) & ag)
        if len(idx) < 2:
            continue
        pos = sites.pos.to_numpy()[idx]
        o = np.argsort(pos)
        p = pos[o]
        near = np.zeros(len(p), bool)
        near[1:] |= np.diff(p) <= window
        near[:-1] |= np.diff(p) <= window
        out[idx[o[near]]] = True
    return out


def _betabinom_logpmf(a, d, al, be):
    return (gammaln(d + 1) - gammaln(a + 1) - gammaln(d - a + 1)
            + betaln(a + al, d - a + be) - betaln(al, be))


def _binom_logpmf(a, d, p):
    return stats.binom.logpmf(a, d, p)


def site_low_vaf_counts(A, D, max_vaf: float = 0.1):
    """Per site: alt and depth summed over cells whose alt fraction is low (likely wild type)."""
    low = (D > 0) & (A <= np.maximum(1, max_vaf * D))
    return (A * low).sum(0), (D * low).sum(0)


def substitution_error_from_sums(subst, a_site, d_site) -> dict:
    rate = (a_site + 0.5) / (d_site + 500)
    good = d_site >= 20
    overall = float(np.median(rate[good])) if good.any() else 1e-3
    out = {}
    subst = pd.Series(np.asarray(subst))
    for t, idx in subst.groupby(subst).groups.items():
        idx = np.asarray(idx)
        ok = d_site[idx] >= 20
        cls = float(np.median(rate[idx][ok])) if ok.any() else overall
        out[t] = min(cls, 3 * overall)
    return out


def substitution_error(sites: pd.DataFrame, A: np.ndarray, D: np.ndarray, max_vaf: float = 0.1) -> dict:
    """Background error per substitution class: the median, over sites of that class, of the alt rate in cells
    whose alt fraction is low (likely wild type). The median keeps a minority of artefact sites (editing,
    hotspots) from inflating the background of their class."""
    out = {}
    low = (D > 0) & (A <= np.maximum(1, max_vaf * D))
    a_site = (A * low).sum(0)
    d_site = (D * low).sum(0)
    rate = (a_site + 0.5) / (d_site + 500)
    good = d_site >= 20
    overall = float(np.median(rate[good])) if good.any() else 1e-3
    for t, idx in sites.groupby("subst").groups.items():
        idx = np.asarray(idx)
        ok = d_site[idx] >= 20
        cls = float(np.median(rate[idx][ok])) if ok.any() else overall
        # a class dominated by one artefact (e.g. T>C full of editing sites) cannot set its own background
        out[t] = min(cls, 3 * overall)
    return out


def call_cells(sites: pd.DataFrame, A: np.ndarray, D: np.ndarray, mut_mean: float = 0.5,
               mut_conc: float = 4.0, prior_strength: float = 2000.0, n_iter: int = 50, min_depth: int = 2,
               learn_vaf: bool = True, vaf_prior: float = 20.0, sub_err: dict | None = None):
    """Per-site EM; returns (posterior matrix cells x sites; NaN = not genotypable (depth < min_depth),
    kept separate from wild type; site table with eps, pi, mutant VAF).

    learn_vaf: the mutant allele fraction of each site is estimated from the cells assigned to the mutant state
    (shrunk toward mut_mean with vaf_prior pseudo-reads), so a cell with few reads is judged against what this
    mutation looks like (e.g. ~0.97 for a homozygous edit, ~0.5 for a heterozygous somatic SNV)."""
    sub_err = sub_err if sub_err is not None else substitution_error(sites, A, D)
    post = np.full(A.shape, np.nan)
    eps_out, pi_out, bg_out, vaf_out = [], [], [], []
    for k in range(A.shape[1]):
        a, d = A[:, k].astype(float), D[:, k].astype(float)
        m = d >= min_depth
        e0 = sub_err.get(sites.subst.iloc[k], 1e-3)
        bg_out.append(e0)
        if m.sum() == 0:
            eps_out.append(e0)
            pi_out.append(0.0)
            vaf_out.append(mut_mean)
            continue
        a, d = a[m], d[m]
        pi, eps, mu = 0.1, e0, mut_mean
        for _ in range(n_iter):
            al, be = mu * mut_conc, (1 - mu) * mut_conc
            l1 = np.log(pi) + _betabinom_logpmf(a, d, al, be)
            l0 = np.log1p(-pi) + _binom_logpmf(a, d, eps)
            mx = np.maximum(l0, l1)
            r = np.exp(l1 - mx) / (np.exp(l0 - mx) + np.exp(l1 - mx))
            pi_new = np.clip((r.sum() + 0.1) / (len(r) + 1.0), 1e-4, 0.999)
            w = 1 - r
            eps_new = ((w * a).sum() + prior_strength * e0) / ((w * d).sum() + prior_strength)
            eps_new = float(np.clip(eps_new, 1e-6, 0.2))
            if learn_vaf:
                mu = float(np.clip(((r * a).sum() + vaf_prior * mut_mean) / ((r * d).sum() + vaf_prior), 0.4, 0.98))
            if abs(pi_new - pi) < 1e-6 and abs(eps_new - eps) < 1e-9:
                pi, eps = pi_new, eps_new
                break
            pi, eps = pi_new, eps_new
        post[m, k] = r
        eps_out.append(eps)
        pi_out.append(pi)
        vaf_out.append(mu)
    s = sites.copy()
    s["error_rate"] = eps_out
    s["background_error"] = bg_out
    s["frac_mutant"] = pi_out
    s["mutant_vaf"] = vaf_out
    return post, s


def genotype_calls(post: np.ndarray, post_call: np.ndarray, threshold: float, frac_mutant,
                   wt_bayes_factor: float = 10.0) -> np.ndarray:
    """Four states per cell and site:
      mutant          callable posterior >= threshold
      wild_type       the cell's own reads favour wild type over the site's mutant model by >= wt_bayes_factor
                      (independent of how common the mutation is at the site)
      uncertain       neither: e.g. 2 reference reads at a site whose mutant VAF is 0.7, or a single alt read
      not_genotypable too few reads
    frac_mutant: per-site mixing proportion from call_cells (to turn posteriors back into likelihood ratios)."""
    pi = np.clip(np.broadcast_to(np.asarray(frac_mutant, float), post.shape), 1e-6, 1 - 1e-6)
    out = np.full(post.shape, "uncertain", dtype=object)
    out[np.isnan(post)] = "not_genotypable"
    fin = np.isfinite(post)
    with np.errstate(divide="ignore", invalid="ignore"):
        bf_wt = pi * (1 - post) / ((1 - pi) * post)       # L(wild type) / L(mutant)
    out[fin & (np.nan_to_num(bf_wt, nan=0, posinf=np.inf) >= wt_bayes_factor)] = "wild_type"
    out[fin & (np.nan_to_num(post_call, nan=0) >= threshold)] = "mutant"
    return out


def callable_posterior(post: np.ndarray, A: np.ndarray, min_cell_alt: int = 2) -> np.ndarray:
    """Posterior used for calls: a cell with fewer than min_cell_alt alt reads cannot be called mutant
    (one read is as likely to be a stray read from another cell - index hopping, plate cross-contamination -
    as a mutant allele). The cell stays genotyped (not NaN): it is called wild type or left uncalled."""
    out = post.copy()
    low = (A < min_cell_alt) & np.isfinite(post)
    out[low] = np.minimum(out[low], 0.0)
    return out


def calibrate_threshold(post: np.ndarray, negative: np.ndarray, usable: np.ndarray, target_fp: float = 0.001,
                        default: float = 0.9):
    """Lowest posterior threshold whose false-call rate in negative-control cells (calls / covered
    cell-site pairs) stays <= target_fp. Returns (threshold, table of threshold vs negative FP rate)."""
    grid = np.array([0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999])
    neg = post[negative][:, usable]
    covered = np.isfinite(neg)
    n = covered.sum()
    rows = []
    for t in grid:
        calls = int((np.nan_to_num(neg, nan=0) >= t).sum())
        rows.append({"threshold": t, "negative_calls": calls, "negative_pairs": int(n),
                     "negative_fp_rate": calls / n if n else np.nan})
    tab = pd.DataFrame(rows)
    ok = tab[tab.negative_fp_rate <= target_fp]
    thr = float(ok.threshold.min()) if len(ok) else float(grid.max())
    return (thr if n else default), tab


def split_check(A: np.ndarray, D: np.ndarray, sites: pd.DataFrame, usable: np.ndarray, threshold: float,
                seed: int = 0):
    """Reproducibility: split each cell's reads at each site into two random halves, call each half
    separately, and report how often a call in one half is supported (alt reads > 0) in the other."""
    rng = np.random.default_rng(seed)
    A1 = rng.hypergeometric(np.maximum(A, 0), np.maximum(D - A, 0), D // 2) if A.size else A
    D1 = D // 2
    A2, D2 = A - A1, D - D1
    p1, _ = call_cells(sites, A1, D1)
    p2, _ = call_cells(sites, A2, D2)
    c1 = (np.nan_to_num(p1, nan=0) >= threshold)[:, usable]
    c2 = (np.nan_to_num(p2, nan=0) >= threshold)[:, usable]
    both = ((D1 > 0) & (D2 > 0))[:, usable]
    supported = ((c1 & both & (A2[:, usable] > 0)).sum() + (c2 & both & (A1[:, usable] > 0)).sum())
    calls = (c1 & both).sum() + (c2 & both).sum()
    return {"half_calls": int(calls), "supported_in_other_half": int(supported),
            "support_rate": float(supported / calls) if calls else np.nan}
