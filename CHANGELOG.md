# Changelog

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
