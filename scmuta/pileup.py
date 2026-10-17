"""Per-cell allele counts with read-level artefact features.

Two passes:
  1. discover: scan every cell, keep positions where a non-reference base passes the read filters
     (pooled over cells) -> candidate sites
  2. count: at candidate sites only, count every base per cell, and for each site sum the read
     features that artefact filters need (strand, distance to read end, distance to a splice
     junction in the read, soft-clip proximity, base quality)

Cells come either from one BAM per cell (plate / Smart-seq2) or from a cell-barcode tag (droplet).
"""
from __future__ import annotations

import os
from collections import defaultdict
from dataclasses import dataclass

import numpy as np
import pandas as pd
import pysam

from .util import log

BASES = "ACGT"
B2I = {b: i for i, b in enumerate(BASES)}


@dataclass
class ReadFilter:
    min_mapq: int = 20
    min_baseq: int = 20
    skip_duplicates: bool = True
    skip_secondary: bool = True
    max_mismatch_frac: float = 0.05     # NM / aligned length; drops mis-mapped reads (paralogs)

    def read_ok(self, r: pysam.AlignedSegment) -> bool:
        if r.is_unmapped or r.mapping_quality < self.min_mapq:
            return False
        if self.skip_secondary and (r.is_secondary or r.is_supplementary):
            return False
        if self.skip_duplicates and r.is_duplicate:
            return False
        if r.has_tag("NM"):
            alen = r.query_alignment_length or 1
            if r.get_tag("NM") / alen > self.max_mismatch_frac:
                return False
        return True


def _junction_and_clip_distance(r: pysam.AlignedSegment):
    """Reference positions of splice junction edges (CIGAR N) and soft-clip edges in this read."""
    edges, clips = [], []
    pos = r.reference_start
    ct = r.cigartuples or []
    for k, (op, ln) in enumerate(ct):
        if op == 3:                      # N: intron
            edges += [pos, pos + ln]
            pos += ln
        elif op in (0, 2, 7, 8):         # M, D, =, X consume reference
            pos += ln
        elif op == 4:                    # S: soft clip at a read end
            clips.append(r.reference_start if k == 0 else pos)
    return edges, clips


def _cells_from_inputs(bams: dict | str, tag: str | None):
    """bams: dict cell -> bam path (plate), or a single BAM path with a barcode tag (droplet)."""
    if isinstance(bams, dict):
        return [(c, p, None) for c, p in bams.items()]
    if not tag:
        raise ValueError("a single BAM needs a cell-barcode tag (e.g. CB)")
    return [(None, bams, tag)]


def _iter_reads(path, region, rf, tag, keep_cells):
    with pysam.AlignmentFile(path) as bam:
        chrom, start, end = region
        for r in bam.fetch(chrom, start, end):
            if not rf.read_ok(r):
                continue
            if tag:
                if not r.has_tag(tag):
                    continue
                cb = r.get_tag(tag)
                if keep_cells is not None and cb not in keep_cells:
                    continue
                yield cb, r
            else:
                yield None, r


def discover(bams, fasta: str, regions, rf: ReadFilter | None = None, tag: str | None = None,
             keep_cells=None, min_alt_reads: int = 2, min_alt_cells: int = 1) -> pd.DataFrame:
    """Pass 1: candidate positions with >= min_alt_reads non-reference reads pooled over cells, seen in
    >= min_alt_cells cells. regions: iterable of (chrom, start, end), 0-based half-open."""
    rf = rf or ReadFilter()
    ref = pysam.FastaFile(fasta)
    alt_reads = defaultdict(int)
    alt_cells = defaultdict(set)
    inputs = _cells_from_inputs(bams, tag)
    for region in regions:
        chrom, start, end = region
        refseq = ref.fetch(chrom, start, end).upper()
        for cell, path, t in inputs:
            for cb, r in _iter_reads(path, region, rf, t, keep_cells):
                cid = cb if t else cell
                q = r.query_qualities
                seq = r.query_sequence
                for qp, rp in r.get_aligned_pairs(matches_only=True):
                    if rp < start or rp >= end or q[qp] < rf.min_baseq:
                        continue
                    b = seq[qp]
                    rb = refseq[rp - start]
                    if b != rb and b in B2I and rb in B2I:
                        key = (chrom, rp, b)
                        alt_reads[key] += 1
                        alt_cells[key].add(cid)
    rows = [(c, p, ref.fetch(c, p, p + 1).upper(), b, n, len(alt_cells[(c, p, b)]))
            for (c, p, b), n in alt_reads.items()
            if n >= min_alt_reads and len(alt_cells[(c, p, b)]) >= min_alt_cells]
    log(f"discover: {len(rows)} candidate alleles from {len(alt_reads)} observed")
    return pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt", "alt_reads", "alt_cells"])


def count(bams, fasta: str, sites: pd.DataFrame, rf: ReadFilter | None = None, tag: str | None = None,
          keep_cells=None, end_window: int = 10, junction_window: int = 10, umi_tag: str | None = None,
          dedup: str = "none"):
    """Pass 2: per-cell ref/alt counts at candidate sites, plus per-site artefact features of alt reads.

    dedup: "none"  - count reads
           "umi"   - count distinct UMIs (umi_tag, e.g. UB) per cell, site and allele (droplet data)
           "start" - count distinct (start, strand) per cell, site and allele; a proxy for molecules in
                     plate data without UMIs (removes PCR duplicates that escaped duplicate marking)
    Returns (cells x sites alt matrix, cells x sites depth matrix (ref+alt), site feature table, cell list).
    """
    if dedup == "umi" and not umi_tag:
        raise ValueError("dedup='umi' needs umi_tag (e.g. UB)")
    seen = set()
    rf = rf or ReadFilter()
    sites = sites.reset_index(drop=True)
    by_pos = defaultdict(list)
    for i, s in sites.iterrows():
        by_pos[(s.chrom, int(s.pos))].append((i, s.ref, s.alt))
    alt = defaultdict(int)      # (cell, site) -> alt reads
    refc = defaultdict(int)     # (cell, site) -> ref reads
    feat = np.zeros((len(sites), 8))  # alt_fwd, alt_rev, ref_fwd, ref_rev, alt_near_end, alt_near_junction,
    #                                   alt_near_clip, alt_baseq_sum
    inputs = _cells_from_inputs(bams, tag)
    # group candidate positions into windows to fetch once per window
    pos_df = sites[["chrom", "pos"]].drop_duplicates().sort_values(["chrom", "pos"])
    windows = []
    for chrom, g in pos_df.groupby("chrom"):
        p = g.pos.to_numpy()
        s0 = p[0]
        prev = p[0]
        for x in p[1:]:
            if x - prev > 5000:
                windows.append((chrom, int(s0), int(prev) + 1))
                s0 = x
            prev = x
        windows.append((chrom, int(s0), int(prev) + 1))
    for cell, path, t in inputs:
        for region in windows:
            for cb, r in _iter_reads(path, region, rf, t, keep_cells):
                cid = cb if t else cell
                q = r.query_qualities
                seq = r.query_sequence
                edges, clips = _junction_and_clip_distance(r)
                qlen = r.query_length
                for qp, rp in r.get_aligned_pairs(matches_only=True):
                    hits = by_pos.get((region[0], rp))
                    if not hits or q[qp] < rf.min_baseq:
                        continue
                    b = seq[qp]
                    near_end = min(qp, qlen - 1 - qp) < end_window
                    near_j = any(abs(rp - e) < junction_window for e in edges)
                    near_c = any(abs(rp - e) < end_window for e in clips)
                    for i, rb, ab in hits:
                        if dedup != "none":
                            mol = (r.get_tag(umi_tag) if r.has_tag(umi_tag) else r.query_name) if dedup == "umi" \
                                else (r.reference_start, r.is_reverse)
                            key = (cid, i, b, mol)
                            if key in seen:
                                continue
                            seen.add(key)
                        if b == ab:
                            alt[(cid, i)] += 1
                            feat[i, 1 if r.is_reverse else 0] += 1
                            feat[i, 4] += near_end
                            feat[i, 5] += near_j
                            feat[i, 6] += near_c
                            feat[i, 7] += q[qp]
                        elif b == rb:
                            refc[(cid, i)] += 1
                            feat[i, 3 if r.is_reverse else 2] += 1
    cells = sorted({k[0] for k in list(alt) + list(refc)})
    ci = {c: j for j, c in enumerate(cells)}
    A = np.zeros((len(cells), len(sites)), np.int32)
    D = np.zeros((len(cells), len(sites)), np.int32)
    for (c, i), n in alt.items():
        A[ci[c], i] += n
        D[ci[c], i] += n
    for (c, i), n in refc.items():
        D[ci[c], i] += n
    ref = pysam.FastaFile(fasta)
    ctx, rep = [], []
    for c_, p_, a_ in zip(sites.chrom, sites.pos, sites.alt):
        left = ref.fetch(c_, max(0, int(p_) - 4), int(p_)).upper()
        right = ref.fetch(c_, int(p_) + 1, int(p_) + 5).upper()
        ctx.append(f"{left}[{ref.fetch(c_, int(p_), int(p_) + 1).upper()}]{right}")
        # alt base extends an adjacent homopolymer run (e.g. TT[C]TT > T): a known sequencing/library artefact
        rep.append(int(left[-2:] == a_ * 2 or right[:2] == a_ * 2))
    f = pd.DataFrame(feat, columns=["alt_fwd", "alt_rev", "ref_fwd", "ref_rev", "alt_near_end",
                                    "alt_near_junction", "alt_near_clip", "alt_baseq_sum"])
    f["context"] = ctx
    f["repeat_extension"] = rep
    log(f"count: {len(cells)} cells x {len(sites)} sites")
    return A, D, pd.concat([sites, f], axis=1), cells


def plate_bams(table: str) -> dict:
    """Read a two-column table (cell, bam) into {cell: bam}."""
    sep = "\t" if table.endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(table, sep=sep)
    return dict(zip(df.iloc[:, 0].astype(str), df.iloc[:, 1].map(os.path.expanduser)))
