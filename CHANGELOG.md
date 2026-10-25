# Changelog

## 0.2.0 (2026-10-25)
- `count()` stores allele counts sparsely (`sparse_out=True`) and fetches reads in windows of at most 200 kb or 64
  candidate sites (`max_window`, `max_sites`), which bounds memory in genes with dense runs of adjacent candidates
  (chr11 of the motor cortex data previously ran out of memory). Counts do not depend on the window size (tested).
- `examples/snrna_pilot.py`: molecule-level germline test (`germline_by_molecules`, `--germline-min-frac`) that
  catches inherited heterozygous variants seen as about one molecule per nucleus.
- `scmuta annotate`: gene spans and regions use all transcripts of a gene; positions outside the chosen
  transcript are reported as `intron` or `exon_other_isoform` instead of intergenic.
- Motor cortex scripts, renamed from the `als_*` drafts: `snrna_summary.py` (detection limit),
  `snrna_candidates.py`, `snrna_site_check.py` (read-level check, optional plots, `--format pdf`) and
  `snrna_classify.py` (final calls, Supplementary Table S1); `mem_probe.py`.
- Fix: the UMI dedup set in `count()` is cleared per site window instead of growing over a whole chromosome,
  which exhausted memory on large droplet data (42k nuclei x 180k sites). Counts are unchanged.
- Head-to-head benchmark against SComatic and scAllele; scoring restricted to a panel of normals or to
  chosen chromosomes; retry on transient read failures.
- `scmuta annotate` and `scmuta call --gtf`: gene, region, consequence, codon and protein change from a GTF
  (checked codon by codon against a brute-force translation on both strands, including split codons).
- Per-cell genotype EM vectorised over all covered cell-site pairs (identical results, ~5x faster on large
  droplet data; per-site reference implementation kept for tests).

## 0.1.0 (2026-10-21)
- Position filters (read end, soft clip) compare alt reads with reference reads at the same site; junction
  window 3 nt; homopolymer extension requires a run >= 4.
- Leave-one-out calibration on control cells, counted only at sites that would be reported.
- Fast counting for plate data (htslib base counts per strand); CIGAR-based read positions for droplet data;
  spliced reads spanning two site windows are no longer counted twice.
- `scmuta call` counts large runs in blocks into sparse matrices and uses the chunked pipeline.
- Genome-wide cell-line benchmark, per-protocol scoring, competitor harness (SComatic, scAllele).

## 0.1.0.dev0
- Prototype: per-site EM genotypes, four states, negative-control calibration, artefact flags, FUS P525L test.
