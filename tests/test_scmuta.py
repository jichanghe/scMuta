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


def _brute_protein(fa_seq, cds, strand, pos, alt):
    from scmuta.annotate import CODON, COMP
    segs = cds if strand == "+" else cds[::-1]
    def cds_seq(seq):
        parts = [seq[s:e] for s, e in cds]
        x = "".join(parts)
        return x if strand == "+" else x[::-1].translate(COMP)
    mut = fa_seq[:pos] + alt + fa_seq[pos + 1:]
    a, b = cds_seq(fa_seq), cds_seq(mut)
    pa = "".join(CODON[a[i:i + 3]] for i in range(0, len(a) - 2, 3))
    pb = "".join(CODON[b[i:i + 3]] for i in range(0, len(b) - 2, 3))
    for i, (x, y) in enumerate(zip(pa, pb)):
        if x != y:
            return f"p.{x}{i + 1}{y}"
    return "syn"


def test_annotation_matches_brute_force(tmp_path):
    import random
    import pysam
    from scmuta.annotate import annotate_sites
    rnd = random.Random(1)
    seq = "".join(rnd.choice("ACGT") for _ in range(600))
    (tmp_path / "r.fa").write_text(">chrA\n" + seq + "\n")
    pysam.faidx(str(tmp_path / "r.fa"))
    # + strand gene, CDS split across two exons with a codon spanning the junction; - strand gene likewise
    gtf = [("chrA", "exon", 101, 140, "+", "G1", "T1"), ("chrA", "CDS", 112, 140, "+", "G1", "T1"),
           ("chrA", "exon", 201, 250, "+", "G1", "T1"), ("chrA", "CDS", 201, 231, "+", "G1", "T1"),
           ("chrA", "exon", 351, 380, "-", "G2", "T2"), ("chrA", "CDS", 351, 371, "-", "G2", "T2"),
           ("chrA", "exon", 451, 500, "-", "G2", "T2"), ("chrA", "CDS", 461, 496, "-", "G2", "T2")]
    with open(tmp_path / "g.gtf", "w") as fh:
        for c, f, s, e, st, g, t in gtf:
            fh.write(f'{c}\tt\t{f}\t{s}\t{e}\t.\t{st}\t.\tgene_id "{g}"; transcript_id "{t}"; gene_name "{g}"; '
                     f'gene_type "protein_coding";\n')
    cds = {"G1": [(111, 140), (200, 231)], "G2": [(350, 371), (460, 496)]}
    strand = {"G1": "+", "G2": "-"}
    pos = list(range(111, 140)) + list(range(200, 231)) + list(range(350, 371)) + list(range(460, 496))
    rows = [(p, b) for p in pos for b in "ACGT" if b != seq[p]]
    sites = pd.DataFrame({"chrom": "chrA", "pos": [p for p, _ in rows], "ref": [seq[p] for p, _ in rows],
                          "alt": [b for _, b in rows]})
    ann = annotate_sites(sites, str(tmp_path / "g.gtf"), str(tmp_path / "r.fa"))
    for (p, b), (_, r) in zip(rows, ann.iterrows()):
        expect = _brute_protein(seq, cds[r.gene], strand[r.gene], p, b)
        got = r.protein_change if r.consequence != "synonymous" else "syn"
        assert got == expect, (p, b, r.gene, got, expect)
    other = annotate_sites(pd.DataFrame({"chrom": ["chrA"] * 4, "pos": [50, 105, 141, 300],
                                         "ref": [seq[50], seq[105], seq[141], seq[300]], "alt": ["N"] * 4}),
                           str(tmp_path / "g.gtf"), str(tmp_path / "r.fa"))
    assert other.region.tolist() == ["intergenic", "5UTR", "splice_site", "intergenic"]


def test_umi_dedup_across_windows(tmp_path):
    """Molecules are collapsed per (cell, UMI, site) even though the dedup set is cleared per window."""
    import random
    import pysam
    from scmuta.pileup import count
    rnd = random.Random(0)
    ref = "".join(rnd.choice("ACGT") for _ in range(400000))
    fa = tmp_path / "r.fa"
    fa.write_text(">chrT\n" + "\n".join(ref[i:i + 60] for i in range(0, len(ref), 60)) + "\n")
    pysam.faidx(str(fa))
    positions = [1000, 1200, 300000, 300050]          # two windows, >5000 nt apart
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"LN": len(ref), "SN": "chrT"}]}
    recs, k = [], 0
    for cell in ("c1", "c2"):
        for pos in positions:
            for umi in range(3):                       # 3 molecules
                for _ in range(4):                     # each with 4 PCR copies
                    k += 1
                    a = pysam.AlignedSegment()
                    start = pos - 30
                    seq = list(ref[start:start + 60])
                    seq[30] = "T" if ref[pos] != "T" else "A"
                    a.query_name, a.query_sequence, a.flag = f"r{k}", "".join(seq), 0
                    a.reference_id, a.reference_start, a.mapping_quality = 0, start, 60
                    a.cigar = [(0, 60)]
                    a.query_qualities = pysam.qualitystring_to_array("I" * 60)
                    a.set_tag("NM", 1)
                    a.set_tag("CB", cell)
                    a.set_tag("UB", f"{cell}_{pos}_{umi}")
                    recs.append(a)
    recs.sort(key=lambda r: r.reference_start)
    bam = tmp_path / "m.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for r in recs:
            out.write(r)
    pysam.index(str(bam))
    sites = pd.DataFrame({"chrom": "chrT", "pos": positions, "ref": [ref[p] for p in positions],
                          "alt": ["T" if ref[p] != "T" else "A" for p in positions]})
    A, D, s, cells = count(str(bam), str(fa), sites, tag="CB", umi_tag="UB", dedup="umi")
    assert cells == ["c1", "c2"]
    assert (A == 3).all(), A
