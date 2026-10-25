#!/usr/bin/env python
"""Final read-level classification of snRNA candidate sites (Supplementary Table S1 of the paper).

Inputs: the annotated candidate table from snrna_candidates.py (scMuta site statistics, gene annotation and
segmental-duplication identity) and site_summary.tsv from snrna_site_check.py (molecule counts in carrier nuclei,
the donor's other nuclei of the same cell type, the same cell type in other donors, and the donor's other cell
types). They are joined on donor, chromosome and position; a table already joined can be given with --sites. Each site gets site flags, a tier and a final call:

  1 Inherited (germline-like)   non-carrier nuclei of the same donor and type carry the allele (>= 15% of
                                molecules, above the other donors' rate, P < 1e-3), or the donor's other cell
                                types do (>= 15% of >= 20 molecules, P < 1e-3)
  2 Candidate somatic           consistent with a clone, no site flag, < 2% alt in other donors from >= 20
                                molecules, and more than one alt molecule per carrier nucleus on recount
  3 Unresolved                  alt above background but < 15%, or too few molecules to tell
  4 Clone-consistent but artefact or weak

  python examples/snrna_classify.py --candidates snrna_candidates.tsv --site-check site_check_all/site_summary.tsv \
      --out classified.tsv \
      --groups 210=C9orf72_FTLD,214=C9orf72_FTLD,224=C9orf72_FTLD,301=Control,302=Control,303=Control
"""
import argparse
import numpy as np
import pandas as pd


def site_flags(m, cluster_nt=25, dense_nt=1000, dense_n=3, recur_nt=10000, segdup_min=0.99):
    flags = [[] for _ in range(len(m))]
    for i, r in m.iterrows():
        same = m[(m.donor == r.donor) & (m.chrom == r.chrom)]
        near = same[(same.index != i) & ((same.pos - r.pos).abs() <= cluster_nt)]
        if len(near):
            flags[i].append(f"clustered ({len(near) + 1} sites within {cluster_nt} nt, same donor)")
        dense = same[(same.pos - r.pos).abs() <= dense_nt]
        if len(dense) >= dense_n and not len(near):
            flags[i].append(f"dense cluster ({len(dense)} candidates within 1 kb, same donor)")
        if r.segdup_identity >= segdup_min:
            flags[i].append(f"segmental duplication, {r.segdup_identity:.1%} identity")
        o = m[(m.donor != r.donor) & (m.chrom == r.chrom) & ((m.pos - r.pos).abs() <= recur_nt)]
        if len(o):
            flags[i].append("locus recurs in donor " + ",".join(sorted(set(o.donor))))
        if r.negative_cells_called > 0:
            flags[i].append(f"also called in {int(r.negative_cells_called)} nuclei of other donors")
        if abs(r.mutant_vaf - 0.4) < 1e-9:
            flags[i].append("mutant VAF at starting value (0.4)")
        if r.vaf_pooled >= 0.2:
            flags[i].append("pooled VAF >= 0.2: possible germline")
    return flags


def final_call(r, other_types_inherited):
    v = r.verdict
    frac = r.noncarrier_alt_fraction
    if v.startswith("germline") or (other_types_inherited and not v.startswith("consistent")):
        return "1 Inherited (germline-like)", ("" if v.startswith("germline") else "also in the donor's other cell types")
    small_confined = v.startswith("confined") and frac < 0.05
    if v.startswith("consistent") or small_confined:
        weak = [x for x, c in [
            ("artefact flag", r.n_flags > 0),
            ("noisy site: >=2% alt in other donors", r.background_alt_fraction_other_donors >= 0.02),
            ("call not replicated: <=1 alt molecule per carrier on recount", r.carriers_alt_molecules <= r.nuclei),
            ("controls uninformative: <20 molecules in other donors", r.other_donors_same_type_molecules < 20)] if c]
        return ("4 Clone-consistent but artefact or weak", "; ".join(weak)) if weak else ("2 Candidate somatic", "")
    note = ("alt in about the same fraction of every cell type: early mosaic or allelic imbalance"
            if v.startswith(("confined", "unresolved")) else "")
    return "3 Unresolved", note


def join_tables(cand, check):
    """Candidate table (0-based pos) joined with snrna_site_check.py's site_summary.tsv (pos_1based)."""
    check = check.copy()
    check["pos"] = check["pos_1based"] - 1
    keep = ["donor", "chrom", "pos"] + [c for c in check.columns if c not in cand.columns]
    m = cand.merge(check[keep], on=["donor", "chrom", "pos"], how="inner")
    if len(m) < len(cand):
        print(f"note: {len(cand) - len(m)} candidate sites have no site-check row and were skipped")
    return m


def classify(m, groups=None):
    m = m.copy().reset_index(drop=True)
    m["donor"] = m["donor"].astype(str)
    for c in ("gene", "other_genes", "protein_change", "region", "consequence"):
        if c in m.columns:
            m[c] = m[c].fillna("")
    if "segdup_identity" not in m.columns:
        m["segdup_identity"] = 0.0
    m["segdup_identity"] = m["segdup_identity"].fillna(0.0)
    m["pos1"] = m.pos + 1
    if groups:
        m["group"] = m.donor.map(groups)
    fl = site_flags(m)
    m["flags"] = ["; ".join(f) for f in fl]
    m["n_flags"] = [len(f) for f in fl]
    m["tier"] = np.where(m.n_flags > 0, "C", np.where(m.nuclei >= 5, "A", "B"))
    xt = ((m.same_donor_other_types_molecules >= 20) & (m.other_types_alt_fraction >= 0.15)
          & (m.other_types_p_above_background < 1e-3))
    res = [final_call(r, x) for (_, r), x in zip(m.iterrows(), xt)]
    m["final"] = [a for a, _ in res]
    m["note"] = [b for _, b in res]
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", help="annotated candidate table (snrna_candidates.py, pos 0-based)")
    ap.add_argument("--site-check", help="site_summary.tsv from snrna_site_check.py")
    ap.add_argument("--sites", help="instead of the two above: a table that already joins them")
    ap.add_argument("--out", default="classified.tsv")
    ap.add_argument("--groups", default="", help="donor=group,... (underscores become spaces)")
    a = ap.parse_args()
    groups = dict(kv.split("=", 1) for kv in a.groups.split(",") if kv) if a.groups else None
    if groups:
        groups = {k: v.replace("_", " ") for k, v in groups.items()}
    if a.sites:
        sites = pd.read_csv(a.sites, sep="\t", dtype={"donor": str})
    elif a.candidates and a.site_check:
        sites = join_tables(pd.read_csv(a.candidates, sep="\t", dtype={"donor": str}, low_memory=False),
                            pd.read_csv(a.site_check, sep="\t", dtype={"donor": str}))
    else:
        ap.error("give --candidates and --site-check, or --sites")
    m = classify(sites, groups)
    m.to_csv(a.out, sep="\t", index=False)
    print(m.final.value_counts().sort_index().to_string())
    if groups:
        print("\nby group:")
        print(pd.crosstab(m.final, m.group).to_string())
    print("\ntiers:", m.tier.value_counts().sort_index().to_dict())
    kinds = ("clustered", "dense cluster", "segmental duplication", "locus recurs", "also called",
             "mutant VAF at starting value", "pooled VAF")
    print("sites per flag:", {k: int(m["flags"].str.contains(k, regex=False).sum()) for k in kinds})
    inh = m[m.final.str.startswith("1")].noncarrier_alt_fraction
    if len(inh):
        print(f"inherited: non-carrier alt fraction median {inh.median():.3f} "
              f"(IQR {inh.quantile(.25):.3f}-{inh.quantile(.75):.3f})")
    cols = [c for c in ["donor", "chrom", "pos1", "subst", "gene", "cell_types", "carriers_alt_molecules",
                        "same_type_noncarriers_alt_molecules", "same_type_noncarriers_molecules",
                        "other_donors_same_type_alt_molecules", "other_donors_same_type_molecules"] if c in m.columns]
    print("\ncandidate somatic:")
    print(m[m.final.str.startswith("2")].sort_values("nuclei", ascending=False)[cols].to_string(index=False))
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
