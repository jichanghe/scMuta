"""Gene and protein-change annotation of SNVs from a GENCODE/Ensembl GTF and the reference FASTA.

For each site (0-based `pos`, `ref`, `alt`) the overlapping genes are found and, per gene, one transcript is used:
MANE_Select > Ensembl_canonical > longest coding sequence > longest transcript. Output columns:
  gene, gene_id, gene_type, strand, transcript_id, region (CDS / 5UTR / 3UTR / noncoding_exon / splice_site /
  intron / intergenic), consequence (missense / synonymous / stop_gained / stop_lost / start_lost / splice_site /
  utr / noncoding / intronic / intergenic), codon_change (ref>alt codon, changed base upper case), protein_change
  (p.V617F) and other_genes (further overlapping genes).
No external tools are needed; a GTF of a few GB is read once, keeping only the chromosomes that carry sites.
"""
from __future__ import annotations

import gzip
import re
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam

CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, x in enumerate(_b):
    for j, y in enumerate(_b):
        for k, z in enumerate(_b):
            CODON[x + y + z] = _aa[16 * i + 4 * j + k]
THREE = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His",
         "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
         "Y": "Tyr", "V": "Val", "*": "Ter"}
COMP = str.maketrans("ACGTacgt", "TGCAtgca")
_ATTR = re.compile(r'(\S+) "([^"]*)"')


def _attrs(s):
    d, tags = {}, []
    for k, v in _ATTR.findall(s):
        if k == "tag":
            tags.append(v)
        else:
            d[k] = v
    d["tags"] = tags
    return d


def read_gtf(path, chroms=None):
    """Genes and their chosen transcript (exons + CDS), for the given chromosomes only."""
    chroms = set(chroms) if chroms is not None else None
    tx = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 9 or f[2] not in ("exon", "CDS") or (chroms is not None and f[0] not in chroms):
                continue
            a = _attrs(f[8])
            tid = a.get("transcript_id")
            if not tid:
                continue
            t = tx.get(tid)
            if t is None:
                t = tx[tid] = {"chrom": f[0], "strand": f[6], "gene": a.get("gene_name", a.get("gene_id", "")),
                               "gene_id": a.get("gene_id", ""), "gene_type": a.get("gene_type", a.get("gene_biotype", "")),
                               "tags": set(a["tags"]), "exons": [], "cds": []}
            (t["exons"] if f[2] == "exon" else t["cds"]).append((int(f[3]) - 1, int(f[4])))     # 0-based half-open
            t["tags"].update(a["tags"])
    by_gene = defaultdict(list)
    for tid, t in tx.items():
        t["id"] = tid
        t["exons"].sort()
        t["cds"].sort()
        by_gene[t["gene_id"]].append(t)

    def rank(t):
        cds_len = sum(e - s for s, e in t["cds"])
        ex_len = sum(e - s for s, e in t["exons"])
        return ("MANE_Select" in t["tags"], "Ensembl_canonical" in t["tags"], cds_len, ex_len)

    genes = []
    for gid, ts in by_gene.items():
        t = max(ts, key=rank)
        # the gene spans all its transcripts; positions outside the chosen transcript but inside another isoform are
        # reported as intronic (or as an exon of another isoform), not intergenic
        g = {**t, "tx_start": t["exons"][0][0], "tx_end": t["exons"][-1][1],
             "start": min(x["exons"][0][0] for x in ts), "end": max(x["exons"][-1][1] for x in ts),
             "other_exons": sorted({e for x in ts if x is not t for e in x["exons"]})}
        genes.append(g)
    return genes


class Annotator:
    def __init__(self, gtf, fasta, chroms=None):
        self.fa = pysam.FastaFile(fasta)
        self.genes = read_gtf(gtf, chroms)
        self.idx = {}
        for c in {g["chrom"] for g in self.genes}:
            gs = sorted((g for g in self.genes if g["chrom"] == c), key=lambda g: g["start"])
            self.idx[c] = (gs, np.array([g["start"] for g in gs]), np.maximum.accumulate([g["end"] for g in gs]))

    def overlapping(self, chrom, pos):
        if chrom not in self.idx:
            return []
        gs, starts, max_end = self.idx[chrom]
        hi = np.searchsorted(starts, pos, side="right")
        lo = np.searchsorted(max_end, pos, side="right")          # genes before lo all end at or before pos
        return [g for g in gs[lo:hi] if g["start"] <= pos < g["end"]]

    def _one(self, g, pos, ref, alt):
        out = {"gene": g["gene"], "gene_id": g["gene_id"], "gene_type": g["gene_type"], "strand": g["strand"],
               "transcript_id": g["id"], "codon_change": "", "protein_change": ""}
        if not (g.get("tx_start", g["start"]) <= pos < g.get("tx_end", g["end"])):
            if any(s <= pos < e for s, e in g.get("other_exons", ())):
                out["region"], out["consequence"] = "exon_other_isoform", "other_isoform"
            else:
                out["region"], out["consequence"] = "intron", "intronic"
            return out
        in_exon = any(s <= pos < e for s, e in g["exons"])
        if not in_exon:
            # within 2 nt of an exon edge, inside the intron = canonical splice site
            near = any((s - 2 <= pos < s) or (e <= pos < e + 2) for s, e in g["exons"])
            out["region"] = "splice_site" if near else "intron"
            out["consequence"] = "splice_site" if near else "intronic"
            return out
        cds = g["cds"]
        if not cds:
            out["region"], out["consequence"] = "noncoding_exon", "noncoding"
            return out
        if not any(s <= pos < e for s, e in cds):
            five = (pos < cds[0][0]) if g["strand"] == "+" else (pos >= cds[-1][1])
            out["region"] = "5UTR" if five else "3UTR"
            out["consequence"] = "utr"
            return out
        out["region"] = "CDS"
        # coding position along the transcript direction
        segs = cds if g["strand"] == "+" else cds[::-1]
        cpos = 0
        for s, e in segs:
            if s <= pos < e:
                cpos += (pos - s) if g["strand"] == "+" else (e - 1 - pos)
                break
            cpos += e - s
        ci = cpos // 3
        # genome positions of the codon's three bases
        gp = []
        acc = 0
        for s, e in segs:
            n = e - s
            for k in range(3):
                want = 3 * ci + k
                if acc <= want < acc + n:
                    off = want - acc
                    gp.append(s + off if g["strand"] == "+" else e - 1 - off)
            acc += n
        if len(gp) != 3:
            out["consequence"] = "coding_incomplete"
            return out
        ref_codon = "".join(self.fa.fetch(g["chrom"], p, p + 1).upper() for p in gp)
        alt_codon = "".join(alt.upper() if p == pos else b for p, b in zip(gp, ref_codon))
        if g["strand"] == "-":
            ref_codon, alt_codon = ref_codon.translate(COMP), alt_codon.translate(COMP)
        a0, a1 = CODON.get(ref_codon, "X"), CODON.get(alt_codon, "X")
        k = gp.index(pos)
        out["codon_change"] = (ref_codon[:k] + ref_codon[k].lower() + ref_codon[k + 1:] + ">" +
                               alt_codon[:k] + alt_codon[k].lower() + alt_codon[k + 1:])
        out["protein_change"] = f"p.{a0}{ci + 1}{a1}" if a0 != a1 else f"p.{a0}{ci + 1}="
        if a0 == a1:
            out["consequence"] = "synonymous"
        elif a1 == "*":
            out["consequence"] = "stop_gained"
        elif a0 == "*":
            out["consequence"] = "stop_lost"
        elif ci == 0 and a0 == "M":
            out["consequence"] = "start_lost"
        else:
            out["consequence"] = "missense"
        return out

    def annotate(self, sites: pd.DataFrame) -> pd.DataFrame:
        order = {"stop_gained": 0, "start_lost": 1, "stop_lost": 2, "missense": 3, "splice_site": 4, "synonymous": 5,
                 "coding_incomplete": 6, "utr": 7, "other_isoform": 7, "noncoding": 8, "intronic": 9}
        rows = []
        for c, p, r, a in zip(sites.chrom, sites.pos.astype(int), sites.ref, sites.alt):
            hits = [self._one(g, p, r, a) for g in self.overlapping(c, p)]
            if not hits:
                rows.append({"gene": "", "gene_id": "", "gene_type": "", "strand": "", "transcript_id": "",
                             "region": "intergenic", "consequence": "intergenic", "codon_change": "",
                             "protein_change": "", "other_genes": ""})
                continue
            # protein-coding genes and the most severe consequence first
            hits.sort(key=lambda h: (h["gene_type"] != "protein_coding", order.get(h["consequence"], 10)))
            best = dict(hits[0])
            best["other_genes"] = ",".join(sorted({h["gene"] for h in hits[1:]} - {best["gene"]}))
            rows.append(best)
        ann = pd.DataFrame(rows, index=sites.index)
        return pd.concat([sites, ann[[c for c in ann.columns if c not in sites.columns]]], axis=1)


def annotate_sites(sites: pd.DataFrame, gtf: str, fasta: str) -> pd.DataFrame:
    return Annotator(gtf, fasta, chroms=set(sites.chrom)).annotate(sites)
