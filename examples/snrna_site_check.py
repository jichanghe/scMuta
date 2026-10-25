#!/usr/bin/env python
"""Read-level evidence and a pyGenomeTracks figure for candidate somatic sites in 10x snRNA-seq (snrna_pilot output).

For each site (donor, chromosome, 1-based position) this
  1. takes the carrier nuclei (called mutant at the site) from <pilot>/<chrom>/somatic_calls.tsv, and the cell types
     of all nuclei from <pilot>/<chrom>/coverage_per_nucleus.tsv;
  2. reads the STARsolo BAMs and reports, per molecule (CB + UB), the base at the site, plus for every alt read its
     distance to the read end, base quality, mapping quality, strand and the number of OTHER mismatches against the
     reference within 25 nt (several mismatches carried together = misalignment / paralog signature);
  3. counts alt / total molecules in three groups: carrier nuclei; non-carrier nuclei of the same donor and cell
     type; nuclei of the same cell type in the other donors (controls);
  4. writes coverage (molecules per nucleus, mean) for the three groups as bedGraphs, a GENCODE slice and an .ini,
     and runs pyGenomeTracks to draw the region with the site marked.

Outputs in --out: <site>.png, <site>.ini, <site>_alt_reads.tsv, and site_summary.tsv for all sites.

  python examples/snrna_site_check.py --pilot ~/snrna_genome --solo STARsolo \
      --fasta GRCh38.primary_assembly.genome.fa --gtf gencode.v50.primary_assembly.annotation.gtf \
      --site 210:chr9:17797007 --site 303:chr8:31701543 --out site_plots
pyGenomeTracks can live in its own conda env: --pgt "conda run -n pgt pyGenomeTracks".
"""
from __future__ import annotations

import argparse
import os
import shlex
import subprocess
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import pysam
from scipy.stats import binomtest

CARRIER, SAME, OTHER, XTYPE = "carriers", "same_type_noncarriers", "other_donors_same_type", "same_donor_other_types"


def read_ok(r, min_mapq):
    return not (r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail) \
        and r.mapping_quality >= min_mapq and r.has_tag("CB") and r.has_tag("UB")


def base_at(r, pos0):
    """Query index aligned to reference position pos0, or None (deletion / splice gap / not covered)."""
    for q, p in r.get_aligned_pairs(matches_only=True):
        if p == pos0:
            return q
    return None


def other_mismatches(r, fa, chrom, pos0, window=25):
    """Mismatches against the reference on this read within +-window of the site, excluding the site itself."""
    pairs = [(q, p) for q, p in r.get_aligned_pairs(matches_only=True) if abs(p - pos0) <= window and p != pos0]
    if not pairs:
        return 0
    lo, hi = min(p for _, p in pairs), max(p for _, p in pairs)
    ref = fa.fetch(chrom, lo, hi + 1).upper()
    seq = r.query_sequence
    return sum(1 for q, p in pairs if seq[q].upper() != ref[p - lo] and ref[p - lo] != "N")


def coverage(bam, chrom, start, end, barcodes, min_mapq):
    """Molecules (CB, UB) covering each base of [start, end), over the given barcodes."""
    mol = defaultdict(set)
    for r in bam.fetch(chrom, start, end):
        if not read_ok(r, min_mapq) or r.get_tag("CB") not in barcodes:
            continue
        key = (r.get_tag("CB"), r.get_tag("UB"))
        for b0, b1 in r.get_blocks():
            for p in range(max(b0, start), min(b1, end)):
                mol[p].add(key)
    cov = np.zeros(end - start)
    for p, s in mol.items():
        cov[p - start] = len(s)
    return cov


def write_bedgraph(path, chrom, start, values):
    with open(path, "w") as fh:
        i, n = 0, len(values)
        while i < n:
            j = i
            while j + 1 < n and values[j + 1] == values[i]:
                j += 1
            fh.write(f"{chrom}\t{start + i}\t{start + j + 1}\t{values[i]:.4g}\n")
            i = j + 1


def gtf_slices(gtf, regions):
    """One pass over the GTF: for each region, the genes and transcripts overlapping it, with ALL their exons / CDS /
    UTRs, including those outside the window (a transcript sliced to the window loses its exons and is drawn as a
    solid block). GENCODE lists each transcript before its exons, so one pass suffices."""
    import re
    tid_re = re.compile(r'transcript_id "([^"]+)"')
    by_chrom = defaultdict(list)
    for k, (c, s, e) in regions.items():
        by_chrom[c].append((k, s, e))
    out = {k: [] for k in regions}
    keep = {k: set() for k in regions}
    parts = ("exon", "CDS", "UTR", "five_prime_utr", "three_prime_utr", "start_codon", "stop_codon")
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if f[0] not in by_chrom or f[2] not in ("gene", "transcript") + parts:
                continue
            s, e = int(f[3]) - 1, int(f[4])
            m = tid_re.search(f[8]) if f[2] != "gene" else None
            for k, rs, re_ in by_chrom[f[0]]:
                if f[2] == "gene":
                    if s < re_ and e > rs:
                        out[k].append(line)
                elif f[2] == "transcript":
                    if s < re_ and e > rs:
                        keep[k].add(m.group(1))
                        out[k].append(line)
                elif m and m.group(1) in keep[k]:
                    out[k].append(line)
    return out


INI = """[spacer]
height = 0.6

[{c}]
file = {pre}_{c}.bedgraph
file_type = bedgraph
title = {t_c}
height = 2.2
color = #2b6cb0
min_value = 0
max_value = {ymax}
fontsize = 9

[spacer]
height = 0.2

[{s}]
file = {pre}_{s}.bedgraph
file_type = bedgraph
title = {t_s}
height = 2.2
color = #7a8691
min_value = 0
max_value = {ymax}
fontsize = 9

[spacer]
height = 0.2

[{o}]
file = {pre}_{o}.bedgraph
file_type = bedgraph
title = {t_o}
height = 2.2
color = #b4bcc4
min_value = 0
max_value = {ymax}
fontsize = 9

[spacer]
height = 0.4

[genes]
file = {pre}_genes.gtf
title =
height = 1.6
style = UCSC
prefered_name = gene_name
merge_transcripts = true
labels = true
labels_in_margin = true
fontsize = 9
color = #4a5560

[spacer]
height = 0.2

[x-axis]
fontsize = 9
where = bottom

[site]
file = {pre}_site.bed
type = vlines
line_width = 0.8
color = #c0392b
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", required=True, help="folder with one snrna_pilot output folder per chromosome")
    ap.add_argument("--solo", required=True, help="folder with one STARsolo folder per donor (STARsolo/<donor>)")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gtf", required=True, help="uncompressed GENCODE GTF")
    ap.add_argument("--site", action="append", default=[], help="DONOR:CHROM:POS (1-based); repeat")
    ap.add_argument("--sites-file", help="table with donor, chrom, pos (0-based, as written by snrna_candidates.py)")
    ap.add_argument("--no-plot", action="store_true", help="read-level counts and verdicts only (fast, for many sites)")
    ap.add_argument("--flank", type=int, default=3000, help="nt drawn on each side of the site")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--pgt", default="pyGenomeTracks", help='command, e.g. "conda run -n pgt pyGenomeTracks"')
    ap.add_argument("--format", default="png", choices=["png", "pdf", "svg"], help="figure format")
    ap.add_argument("--out", default="site_plots")
    a = ap.parse_args()
    if not a.site and not a.sites_file:
        ap.error("give --site and/or --sites-file")
    os.makedirs(a.out, exist_ok=True)
    pilot = os.path.expanduser(a.pilot)
    fa = pysam.FastaFile(a.fasta)

    sites = []
    for s in a.site:
        d, c, p = s.split(":")
        p1 = int(p.replace(",", ""))
        sites.append((d, c, p1 - 1, f"{d}_{c}_{p1}"))
    if a.sites_file:
        t = pd.read_csv(a.sites_file, sep="\t", dtype={"donor": str})
        for d, c, p0 in zip(t.donor, t.chrom, t.pos.astype(int)):
            sites.append((d, c, p0, f"{d}_{c}_{p0 + 1}"))
    regions = {k: (c, max(0, p0 - a.flank), p0 + a.flank + 1) for d, c, p0, k in sites}
    if not a.no_plot:
        print("reading GTF once for all regions ...", flush=True)
        gtf = gtf_slices(a.gtf, regions)

    meta_cache, bams, rows = {}, {}, []
    for d, chrom, p0, key in sites:
        if chrom not in meta_cache:
            meta = pd.read_csv(os.path.join(pilot, chrom, "coverage_per_nucleus.tsv"), sep="\t", dtype={"donor": str})
            calls = pd.read_csv(os.path.join(pilot, chrom, "somatic_calls.tsv"), sep="\t", dtype={"donor": str},
                                low_memory=False)
            meta_cache[chrom] = (meta, calls)
        meta, calls = meta_cache[chrom]
        here = calls[(calls.donor == d) & (calls.chrom == chrom) & (calls.pos == p0)]
        if here.empty:
            print(f"{key}: no calls at this position in {chrom}/somatic_calls.tsv - skipped (check donor and 1-based position)")
            continue
        ref, alt = here.ref.iloc[0], here.alt.iloc[0]
        carriers = set(here.cell.str.split(":", n=1).str[1])
        ctype = Counter(meta.set_index("cell").loc[here.cell, "cell_type"]).most_common(1)[0][0]
        same = set(meta[(meta.donor == d) & (meta.cell_type == ctype)].barcode) - carriers
        others = {od: set(g.barcode) for od, g in meta[(meta.donor != d) & (meta.cell_type == ctype)].groupby("donor")}
        xtype = set(meta[(meta.donor == d) & (meta.cell_type != ctype) & (meta.cell_type != "unassigned")].barcode) - carriers
        groups = {CARRIER: {d: carriers}, SAME: {d: same}, OTHER: others, XTYPE: {d: xtype}}
        n_nuc = {CARRIER: len(carriers), SAME: len(same), OTHER: sum(len(v) for v in others.values()), XTYPE: len(xtype)}

        c0, c1 = regions[key][1], regions[key][2]
        covs, alleles, alt_reads = {}, {}, []
        for g, per_donor in groups.items():
            cov = np.zeros(c1 - c0)
            mol_base = defaultdict(Counter)
            for dd, bcs in per_donor.items():
                if dd not in bams:
                    bams[dd] = pysam.AlignmentFile(os.path.join(a.solo, dd, "Aligned.sortedByCoord.out.bam"))
                bam = bams[dd]
                if not a.no_plot:
                    cov += coverage(bam, chrom, c0, c1, bcs, a.min_mapq)
                for r in bam.fetch(chrom, p0, p0 + 1):
                    if not read_ok(r, a.min_mapq) or r.get_tag("CB") not in bcs:
                        continue
                    q = base_at(r, p0)
                    if q is None:
                        continue
                    b = r.query_sequence[q].upper()
                    mol_base[(dd, r.get_tag("CB"), r.get_tag("UB"))][b] += 1
                    if b == alt:
                        L = r.query_length
                        alt_reads.append({"group": g, "donor": dd, "barcode": r.get_tag("CB"), "umi": r.get_tag("UB"),
                                          "strand": "-" if r.is_reverse else "+", "mapq": r.mapping_quality,
                                          "base_quality": r.query_qualities[q] if r.query_qualities is not None else None,
                                          "dist_to_read_end": min(q, L - 1 - q),
                                          "other_mismatches_25nt": other_mismatches(r, fa, chrom, p0),
                                          "soft_clipped": "S" in (r.cigarstring or "")})
            covs[g] = cov / max(n_nuc[g], 1)
            mol = {k: v.most_common(1)[0][0] for k, v in mol_base.items()}
            alleles[g] = (sum(1 for b in mol.values() if b == alt), len(mol),
                          len({k[:2] for k, b in mol.items() if b == alt}))

        ar = pd.DataFrame(alt_reads)
        ar.to_csv(os.path.join(a.out, f"{key}_alt_reads.tsv"), sep="\t", index=False)
        car = ar[ar.group == CARRIER] if len(ar) else ar
        row = {"site": key, "donor": d, "chrom": chrom, "pos_1based": p0 + 1, "change": f"{ref}>{alt}", "cell_type": ctype}
        for g in (CARRIER, SAME, OTHER, XTYPE):
            k, n, nn = alleles[g]
            row.update({f"{g}_nuclei": n_nuc[g], f"{g}_alt_molecules": k, f"{g}_molecules": n, f"{g}_alt_nuclei": nn})
        if len(car):
            row.update({"carrier_alt_reads": len(car), "carrier_alt_umis": car[["barcode", "umi"]].drop_duplicates().shape[0],
                        "alt_reads_plus_strand": int((car.strand == "+").sum()),
                        "alt_reads_with_other_mismatch": int((car.other_mismatches_25nt > 0).sum()),
                        "median_dist_to_read_end": float(car.dist_to_read_end.median()),
                        "median_base_quality": float(car.base_quality.median()),
                        "min_mapq": int(car.mapq.min()),
                        "alt_reads_soft_clipped": int(car.soft_clipped.sum())})
        # Germline or somatic? The carrier nuclei were chosen because they show the alt allele, so the test uses the
        # NON-carrier nuclei of the same donor and cell type. Under a small clone they show the alt allele only at the
        # background rate (other donors' nuclei); an inherited heterozygous variant is in every nucleus, so with about one
        # molecule per nucleus a third to a half of their molecules carry it (less than half with allelic imbalance).
        ks, ns = alleles[SAME][0], alleles[SAME][1]
        bg = max(alleles[OTHER][0] / max(alleles[OTHER][1], 1), 0.002)
        frac = ks / max(ns, 1)
        implied = min(1.0, 2 * frac)            # fraction of non-carrier nuclei that would have to carry it
        p_non = binomtest(ks, ns, bg, alternative="greater").pvalue if ns else float("nan")
        if ns < 20:
            verdict = "too few molecules in non-carrier nuclei to tell"
        elif p_non < 1e-3 and frac >= 0.15:
            verdict = "germline-like: non-carrier nuclei carry it too"
        elif p_non < 1e-3:
            verdict = "unresolved: alt above background in non-carriers"
        else:
            verdict = "consistent with a clone"
        # an inherited variant is also in the donor's other cell types that express the gene; if those are well covered
        # and clean, the allele is confined to one lineage: a large clone (or cell-type-specific allelic expression)
        kx, nx = alleles[XTYPE][0], alleles[XTYPE][1]
        if verdict.startswith(("germline", "unresolved")) and nx >= 20 and kx / nx < 0.05 and \
                binomtest(kx, nx, bg, alternative="greater").pvalue >= 1e-3:
            verdict = "confined to one cell type: possible large clone"
        row.update({"noncarrier_alt_fraction": round(frac, 4), "background_alt_fraction_other_donors": round(bg, 4),
                    "p_noncarrier_alt_above_background": p_non, "implied_fraction_of_noncarrier_nuclei_carrying": round(implied, 3),
                    "verdict": verdict,
                    # same donor, other cell types: an inherited variant is in every cell type that expresses the gene;
                    # a clone is usually confined to one lineage
                    "other_types_alt_fraction": round(alleles[XTYPE][0] / max(alleles[XTYPE][1], 1), 4),
                    "other_types_p_above_background": binomtest(alleles[XTYPE][0], alleles[XTYPE][1], bg, alternative="greater").pvalue
                    if alleles[XTYPE][1] else float("nan")})
        rows.append(row)
        if a.no_plot:
            print(f"{key}: {ctype}, non-carrier alt {ks}/{ns} ({frac:.0%}) vs background {bg:.2%} -> {verdict}; "
                  f"same donor, other cell types {alleles[XTYPE][0]}/{alleles[XTYPE][1]}", flush=True)
            continue

        # files for pyGenomeTracks
        pre = os.path.join(a.out, key)
        for g in (CARRIER, SAME, OTHER):
            write_bedgraph(f"{pre}_{g}.bedgraph", chrom, c0, covs[g])
        with open(f"{pre}_genes.gtf", "w") as fh:
            fh.writelines(gtf[key])
        with open(f"{pre}_site.bed", "w") as fh:
            fh.write(f"{chrom}\t{p0}\t{p0 + 1}\t{key}\n")
        ymax = max(float(max(v.max() for v in covs.values())), 0.01) * 1.05
        t = lambda g, lab: f"{lab} (n={n_nuc[g]}): {alleles[g][0]}/{alleles[g][1]} {alt}"
        ini = INI.format(c=CARRIER, s=SAME, o=OTHER, pre=os.path.basename(pre), ymax=f"{ymax:.3g}",
                         t_c=t(CARRIER, f"Carrier {ctype}, donor {d}"),
                         t_s=t(SAME, f"Other {ctype}, donor {d}"),
                         t_o=t(OTHER, f"{ctype}, other 5 donors" if len(others) == 5 else f"{ctype}, other donors"))
        with open(f"{pre}.ini", "w") as fh:
            fh.write(ini)
        cmd = shlex.split(a.pgt) + ["--tracks", f"{key}.ini", "--region", f"{chrom}:{c0 + 1}-{c1}",
                                    "--outFileName", f"{key}.{a.format}", "--dpi", "200", "--width", "24", "--trackLabelFraction", "0.3",
                                    "--title", f"{chrom}:{p0 + 1:,} {ref}>{alt}, donor {d} (molecules per nucleus, mean)",
                                    "--fontSize", "10"]
        res = subprocess.run(cmd, cwd=a.out, capture_output=True, text=True)
        status = "ok" if res.returncode == 0 and os.path.exists(f"{pre}.{a.format}") else "FAILED"
        print(f"{key}: {ctype}, carriers {alleles[CARRIER][0]}/{alleles[CARRIER][1]} alt molecules, "
              f"same-type non-carriers {alleles[SAME][0]}/{alleles[SAME][1]}, other donors {alleles[OTHER][0]}/{alleles[OTHER][1]}; "
              f"plot {status}", flush=True)
        if status != "ok":
            print(res.stderr[-1500:])

    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(a.out, "site_summary.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 250)
    if len(summ) <= 6:
        print("\n" + summ.T.to_string(header=False))
    else:
        cols = ["site", "cell_type", "change", "carriers_nuclei", "carriers_alt_molecules", "carriers_molecules",
                "same_type_noncarriers_alt_molecules", "same_type_noncarriers_molecules", "other_donors_same_type_alt_molecules",
                "other_donors_same_type_molecules", "noncarrier_alt_fraction", "same_donor_other_types_alt_molecules",
                "same_donor_other_types_molecules", "verdict"]
        print("\n" + summ[cols].to_string(index=False))
    if "verdict" in summ:
        print("\nverdicts:\n" + summ.verdict.value_counts().to_string())
    print(f"\nwrote {a.out}/site_summary.tsv, <site>_alt_reads.tsv and <site>.{a.format}")


if __name__ == "__main__":
    main()
