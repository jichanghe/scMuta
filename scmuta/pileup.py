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
                if r.has_tag("NM") and r.get_tag("NM") == 0:
                    continue                      # no mismatch in this read: nothing to discover (most reads)
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
          keep_cells=None, end_window: int = 10, junction_window: int = 3, umi_tag: str | None = None,
          dedup: str = "none", fast: bool = True, sparse_out: bool = False, max_window: int = 200_000, max_sites: int = 64):
    """Pass 2: per-cell ref/alt counts at candidate sites, plus per-site artefact features of alt reads.

    dedup: "none"  - count reads
           "umi"   - count distinct UMIs (umi_tag, e.g. UB) per cell, site and allele (droplet data)
           "start" - count distinct (start, strand) per cell, site and allele; a proxy for molecules in
                     plate data without UMIs (removes PCR duplicates that escaped duplicate marking)
    Returns (cells x sites alt matrix, cells x sites depth matrix (ref+alt), site feature table, cell list).
    max_window caps the length of one fetch window. Windows join candidate sites less than 5 kb apart; on a dense
    chromosome they would otherwise merge into a few windows tens of Mb long, and the per-window duplicate set,
    cleared only between windows, would grow until memory ran out. Each site is counted only in the window that holds
    it, so cutting windows anywhere leaves the counts unchanged. max_sites caps the candidate positions per window for
    the same reason: in runs of adjacent candidates inside very deep genes (MALAT1, NEAT1 in single-nucleus data)
    one read covers dozens of sites and the duplicate set grows with reads x sites; smaller windows re-read those
    reads a few more times but keep memory bounded.
    sparse_out=True returns the two matrices as scipy CSR (int32). Use it for droplet data: a dense cells x sites
    table for ~10^4 nuclei and ~10^5 candidate sites takes tens of GB, although almost every entry is zero.
    Plate input (one BAM per cell) uses a fast path: base counts per strand come from htslib
    (count_coverage), and only reads with a mismatch are walked in Python for the artefact features.
    With dedup="start" the fast path keeps the first read per (start, strand) in each cell (the slow path
    deduplicates per site and allele; the two differ only when duplicates disagree at a site).
    """
    if dedup == "umi" and not umi_tag:
        raise ValueError("dedup='umi' needs umi_tag (e.g. UB)")
    if fast and isinstance(bams, dict) and dedup in ("none", "start"):
        return _count_plate_fast(bams, fasta, sites, rf or ReadFilter(), end_window, junction_window, dedup)
    seen = set()
    rf = rf or ReadFilter()
    sites = sites.reset_index(drop=True)
    by_pos = defaultdict(list)
    for i, s in sites.iterrows():
        by_pos[(s.chrom, int(s.pos))].append((i, s.ref, s.alt))
    alt = defaultdict(int)      # (cell, site) -> alt reads, current window
    refc = defaultdict(int)     # (cell, site) -> ref reads, current window
    abuf, rbuf, cell_ix = [], [], {}   # finished windows as (cell index, site, count) arrays; cell -> first-seen index
    feat = np.zeros((len(sites), 11))  # alt_fwd, alt_rev, ref_fwd, ref_rev, alt_near_end, alt_near_junction,
    #                                    alt_near_clip, alt_baseq_sum, ref_near_end, ref_near_junction, ref_near_clip
    inputs = _cells_from_inputs(bams, tag)
    # group candidate positions into windows to fetch once per window
    pos_df = sites[["chrom", "pos"]].drop_duplicates().sort_values(["chrom", "pos"])
    windows = []
    for chrom, g in pos_df.groupby("chrom"):
        p = g.pos.to_numpy()
        s0 = p[0]
        prev = p[0]
        n_in = 1
        for x in p[1:]:
            if x - prev > 5000 or x - s0 >= max_window or n_in >= max_sites:
                windows.append((chrom, int(s0), int(prev) + 1))
                s0 = x
                n_in = 0
            prev = x
            n_in += 1
        windows.append((chrom, int(s0), int(prev) + 1))
    sorted_pos = {c: np.array(sorted({p for (cc, p) in by_pos if cc == c})) for c in sites.chrom.unique()}
    for cell, path, t in inputs:
        for region in windows:
            sp = sorted_pos[region[0]]
            # each site is counted only in the window holding it, so the dedup set (keyed by site) can be cleared
            # per window; windows are capped at max_window, otherwise on dense chromosomes the set would grow to one
            # entry per (cell, UMI, site, allele) over tens of Mb and exhaust memory on large droplet data
            seen.clear()
            for cb, r in _iter_reads(path, region, rf, t, keep_cells):
                # candidate sites inside this read's aligned blocks (most reads in a window cover none): walk the
                # CIGAR to get their read positions directly; sites outside this window are counted in their own
                # window (a spliced read can reach the next one), so they are skipped here
                pairs = []
                rpos, qpos = r.reference_start, 0
                for op, ln in r.cigartuples or ():
                    if op in (0, 7, 8):
                        lo, hi = max(rpos, region[1]), min(rpos + ln, region[2])
                        if hi > lo:
                            j0, j1 = np.searchsorted(sp, lo), np.searchsorted(sp, hi)
                            for x in sp[j0:j1]:
                                pairs.append((qpos + int(x) - rpos, int(x)))
                        rpos += ln
                        qpos += ln
                    elif op in (1, 4):
                        qpos += ln
                    elif op in (2, 3):
                        rpos += ln
                if not pairs:
                    continue
                cid = cb if t else cell
                q = r.query_qualities
                seq = r.query_sequence
                edges, clips = _junction_and_clip_distance(r)
                qlen = r.query_length
                for qp, rp in pairs:
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
                            feat[i, 8] += near_end
                            feat[i, 9] += near_j
                            feat[i, 10] += near_c
            # move this window's counts into compact arrays and empty the dicts: a site belongs to exactly one
            # window, so nothing counted later can add to these entries (a Python dict entry costs ~200 bytes,
            # an array entry 12; on a dense chromosome with ~10^4 nuclei the difference is tens of GB)
            for d_, buf in ((alt, abuf), (refc, rbuf)):
                if d_:
                    n = len(d_)
                    buf.append((np.fromiter((cell_ix.setdefault(c, len(cell_ix)) for c, _ in d_), np.int64, n),
                                np.fromiter((i for _, i in d_), np.int64, n),
                                np.fromiter(d_.values(), np.int32, n)))
                    d_.clear()
    cells = sorted(cell_ix)
    order = np.empty(len(cell_ix), np.int64)              # first-seen index -> position in sorted cell list
    for j, c in enumerate(cells):
        order[cell_ix[c]] = j

    def _cat(buf):
        if not buf:
            return np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0, np.int32)
        return order[np.concatenate([b[0] for b in buf])], np.concatenate([b[1] for b in buf]), \
            np.concatenate([b[2] for b in buf])
    ar, ac, av = _cat(abuf)
    rr, rc, rv = _cat(rbuf)
    del abuf, rbuf
    shape = (len(cells), len(sites))
    if sparse_out:
        from scipy import sparse
        A = sparse.csr_matrix((av, (ar, ac)), shape=shape, dtype=np.int32)
        D = (A + sparse.csr_matrix((rv, (rr, rc)), shape=shape, dtype=np.int32)).tocsr()
    else:
        A = np.zeros(shape, np.int32)
        np.add.at(A, (ar, ac), av)
        D = A.copy()
        np.add.at(D, (rr, rc), rv)
    ref = pysam.FastaFile(fasta)
    ctx, rep = [], []
    for c_, p_, a_ in zip(sites.chrom, sites.pos, sites.alt):
        left = ref.fetch(c_, max(0, int(p_) - 4), int(p_)).upper()
        right = ref.fetch(c_, int(p_) + 1, int(p_) + 5).upper()
        ctx.append(f"{left}[{ref.fetch(c_, int(p_), int(p_) + 1).upper()}]{right}")
        # alt base extends an adjacent homopolymer run (e.g. TT[C]TT > T): a known sequencing/library artefact
        rep.append(int(left[-2:] == a_ * 2 or right[:2] == a_ * 2))
    f = pd.DataFrame(feat, columns=FEATURES)
    f["context"] = ctx
    f["repeat_extension"] = rep
    log(f"count: {len(cells)} cells x {len(sites)} sites")
    return A, D, pd.concat([sites, f], axis=1), cells


FEATURES = ["alt_fwd", "alt_rev", "ref_fwd", "ref_rev", "alt_near_end", "alt_near_junction", "alt_near_clip",
            "alt_baseq_sum", "ref_near_end", "ref_near_junction", "ref_near_clip"]


def _windows(sites, gap=5000, max_span=1_000_000):
    pos_df = sites[["chrom", "pos"]].drop_duplicates().sort_values(["chrom", "pos"])
    out = []
    for chrom, g in pos_df.groupby("chrom", sort=False):
        p = g.pos.to_numpy()
        s0 = prev = p[0]
        for x in p[1:]:
            if x - prev > gap or x - s0 > max_span:
                out.append((chrom, int(s0), int(prev) + 1))
                s0 = x
            prev = x
        out.append((chrom, int(s0), int(prev) + 1))
    return out


def _count_plate_fast(bams, fasta, sites, rf, end_window, junction_window, dedup):
    sites = sites.reset_index(drop=True)
    n = len(sites)
    chrom_a = sites.chrom.to_numpy()
    pos_a = sites.pos.to_numpy().astype(np.int64)
    ref_i = np.array([B2I.get(str(b).upper(), -1) for b in sites.ref])
    alt_i = np.array([B2I.get(str(b).upper(), -1) for b in sites.alt])
    windows = _windows(sites)
    win_sites = []
    for chrom, ws, we in windows:
        idx = np.flatnonzero((chrom_a == chrom) & (pos_a >= ws) & (pos_a < we))
        win_sites.append(idx)
    cell_list = list(bams)
    A = np.zeros((len(cell_list), n), np.int32)
    R = np.zeros((len(cell_list), n), np.int32)
    feat = np.zeros((n, 11))
    for ci, cell in enumerate(cell_list):
        with pysam.AlignmentFile(bams[cell]) as bam:
            for (chrom, ws, we), idx in zip(windows, win_sites):
                if len(idx) == 0:
                    continue
                off = pos_a[idx] - ws
                ok = (ref_i[idx] >= 0) & (alt_i[idx] >= 0) & (ref_i[idx] != alt_i[idx])
                for rev in (False, True):
                    seen = set()

                    def keep(r, rev=rev, seen=seen):
                        if r.is_reverse != rev or not rf.read_ok(r):
                            return False
                        if dedup == "start":
                            k = r.reference_start
                            if k in seen:
                                return False
                            seen.add(k)
                        return True
                    cov = np.asarray(bam.count_coverage(chrom, ws, we, quality_threshold=rf.min_baseq,
                                                      read_callback=keep), dtype=np.int64)
                    a = np.where(ok, cov[np.clip(alt_i[idx], 0, 3), off], 0).astype(np.int32)
                    r_ = np.where(ok, cov[np.clip(ref_i[idx], 0, 3), off], 0).astype(np.int32)
                    A[ci, idx] += a
                    R[ci, idx] += r_
                    feat[idx, 1 if rev else 0] += a
                    feat[idx, 3 if rev else 2] += r_
                # read-level features of alt AND ref bases (alt reads near read ends / junctions are only suspicious
                # when the reference reads at the same site are not): walk the CIGAR, not every aligned pair
                sp = pos_a[idx]
                order = np.argsort(sp)
                sp_sorted = sp[order]
                idx_sorted = idx[order]
                seen = set()
                for r in bam.fetch(chrom, ws, we):
                    if not rf.read_ok(r):
                        continue
                    if dedup == "start":
                        k = (r.reference_start, r.is_reverse)
                        if k in seen:
                            continue
                        seen.add(k)
                    hits = []
                    rpos, qpos = r.reference_start, 0
                    for op, ln in r.cigartuples or ():
                        if op in (0, 7, 8):
                            j0, j1 = np.searchsorted(sp_sorted, rpos), np.searchsorted(sp_sorted, rpos + ln)
                            for j in range(j0, j1):
                                hits.append((int(sp_sorted[j]), qpos + int(sp_sorted[j]) - rpos, int(idx_sorted[j])))
                            rpos += ln
                            qpos += ln
                        elif op in (1, 4):
                            qpos += ln
                        elif op in (2, 3):
                            rpos += ln
                    if not hits:
                        continue
                    q = r.query_qualities
                    seq = r.query_sequence
                    qlen = r.query_length
                    edges, clips = _junction_and_clip_distance(r)
                    for rp, qp, i in hits:
                        if q[qp] < rf.min_baseq:
                            continue
                        b = B2I.get(seq[qp], -1)
                        if b < 0:
                            continue
                        if b == alt_i[i]:
                            o = 4
                        elif b == ref_i[i]:
                            o = 8
                        else:
                            continue
                        feat[i, o] += min(qp, qlen - 1 - qp) < end_window
                        feat[i, o + 1] += any(abs(rp - e) < junction_window for e in edges)
                        feat[i, o + 2] += any(abs(rp - e) < end_window for e in clips)
                        if o == 4:
                            feat[i, 7] += q[qp]
    D = A + R
    keep_rows = np.flatnonzero(D.sum(1) > 0)
    cells = sorted(cell_list[k] for k in keep_rows)
    rows = [cell_list.index(c) for c in cells]
    A, D = A[rows], D[rows]
    ref = pysam.FastaFile(fasta)
    ctx, rep = [], []
    for c_, p_, a_ in zip(sites.chrom, sites.pos, sites.alt):
        left = ref.fetch(c_, max(0, int(p_) - 4), int(p_)).upper()
        right = ref.fetch(c_, int(p_) + 1, int(p_) + 5).upper()
        ctx.append(f"{left}[{ref.fetch(c_, int(p_), int(p_) + 1).upper()}]{right}")
        rep.append(int(left[-2:] == a_ * 2 or right[:2] == a_ * 2))
    f = pd.DataFrame(feat, columns=FEATURES)
    f["context"] = ctx
    f["repeat_extension"] = rep
    log(f"count: {len(cells)} cells x {n} sites")
    return A, D, pd.concat([sites, f], axis=1), cells


def plate_bams(table: str) -> dict:
    """Read a two-column table (cell, bam) into {cell: bam}."""
    sep = "\t" if table.endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(table, sep=sep)
    return dict(zip(df.iloc[:, 0].astype(str), df.iloc[:, 1].map(os.path.expanduser)))
