"""Tiny reference + one BAM per cell with a known somatic SNV, an artefact near a splice junction,
and an editing-like site present in every cell, for integration tests."""
import os
import random

import pysam


def build(d, n_mut=8, n_wt=8, n_ctrl=8, reads=12, seed=0):
    rnd = random.Random(seed)
    os.makedirs(d, exist_ok=True)
    ref = "".join(rnd.choice("ACGT") for _ in range(3000))
    ref = ref[:500] + "A" + ref[501:]           # site 500: somatic A>T
    ref = ref[:1500] + "A" + ref[1501:]         # site 1500: editing-like A>G in every cell
    fa = os.path.join(d, "ref.fa")
    with open(fa, "w") as fh:
        fh.write(">chrT\n" + "\n".join(ref[i:i + 60] for i in range(0, len(ref), 60)) + "\n")
    pysam.faidx(fa)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"LN": len(ref), "SN": "chrT"}]}
    cells = {}
    groups = [("mut", n_mut), ("wt", n_wt), ("ctrl", n_ctrl)]
    for g, n in groups:
        for i in range(n):
            name = f"{g}{i}"
            path = os.path.join(d, f"{name}.bam")
            recs = []

            def add(start, seq, cigar, rev=False, k=[0]):
                k[0] += 1
                a = pysam.AlignedSegment()
                a.query_name = f"{name}_r{k[0]}"
                a.query_sequence = seq
                a.flag = 16 if rev else 0
                a.reference_id = 0
                a.reference_start = start
                a.mapping_quality = 60
                a.cigar = cigar
                a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
                a.set_tag("NM", 0)
                recs.append(a)

            for j in range(reads):
                s = rnd.randint(440, 490)
                seq = list(ref[s:s + 80])
                if g == "mut" and j % 4 in (0, 1):
                    seq[500 - s] = "T"
                add(s, "".join(seq), [(0, 80)], rev=j % 2 == 1)
            for j in range(reads):      # spliced reads with a mismatch right at the junction edge (artefact)
                s = 960
                seq = list(ref[s:s + 40] + ref[1100:1140])
                if j % 3 == 0:
                    seq[39] = "C" if seq[39] != "C" else "G"
                add(s, "".join(seq), [(0, 40), (3, 100), (0, 40)], rev=j % 2 == 1)
            for j in range(reads):
                s = rnd.randint(1450, 1490)
                seq = list(ref[s:s + 80])
                if j % 6 in (0, 1) and (i + j) % 2 == 0:
                    seq[1500 - s] = "G"
                add(s, "".join(seq), [(0, 80)], rev=j % 2 == 1)
            recs.sort(key=lambda r: r.reference_start)
            with pysam.AlignmentFile(path, "wb", header=header) as out:
                for r in recs:
                    out.write(r)
            pysam.index(path)
            cells[name] = path
    with open(os.path.join(d, "cells.tsv"), "w") as fh:
        fh.write("cell\tbam\n" + "".join(f"{c}\t{p}\n" for c, p in cells.items()))
    with open(os.path.join(d, "negative.txt"), "w") as fh:
        fh.write("".join(f"ctrl{i}\n" for i in range(n_ctrl)))
    with open(os.path.join(d, "regions.bed"), "w") as fh:
        fh.write("chrT\t0\t3000\n")
    return fa, cells
