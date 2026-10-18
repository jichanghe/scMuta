# scMuta (working name)

**Per-cell SNV calling from single-cell sequencing that tells you how many calls are false.**
Status: **0.1.0.dev0, a prototype.** Tested on simulations and a small BAM fixture only; not yet on real data.

Single-cell variant callers are rarely checked for calibration on real cells: per-cell calls are dominated
by sequencing errors, RNA editing, alignment artefacts and allelic dropout, and most tools report hard calls
without a measured false-call rate. scMuta follows the scSplicer approach: combine the best existing idea
for each step, then measure false calls directly on real cells that cannot carry the variants.

## What it does

| Step | Idea | Inspired by |
|---|---|---|
| Candidates | Pooled pileup over cells; reads filtered by MAPQ, base quality, duplicates, mismatch load | cellsnp-lite, pooled callers |
| Molecules, not reads | `--dedup umi` (droplet) or `start` (plate) | mtDNA single-molecule critique |
| Artefact flags | strand bias, read-end position, **splice-junction proximity**, soft-clip, base quality, repeat extension (YY[X]YY > Y), known RNA-editing sites, germline-like sites | SComatic, Wiens et al. 2024, scSplicer |
| Per-cell genotype | Mixture per site: wild type = Binomial(error), mutant = Beta-Binomial (allelic imbalance / dropout); error shrunk to a robust per-substitution background; EM | SComatic, scMitoMut, SC-BIG |
| Four states | mutant / wild type (reads favour wild type >= 10:1) / **uncertain** / **not genotypable**; mutant VAF learned per site | gap in current tools |
| **Negative-control calibration** | Cells that cannot carry the variants (another line, donor, unedited cells) screen out sites with excess alt reads (editing, hotspots, germline) and set the posterior threshold; cross-fitted so the reported false-call rate is measured on cells that did not screen those sites | SCEPTRE-style empirical nulls; scSplicer |
| Reproducibility | Each cell's reads split in two; calls in one half checked for support in the other | scSplicer split checks |

## Simulation benchmark

`python examples/benchmark_sim.py` (150 case + 100 negative-control cells; 60 somatic sites in 3 clones;
editing, germline, strand, junction, hotspot and noise sites; depth with many zeros; 30% allelic dropout;
mean of 5 seeds):

| Method | Cell-level precision | Recall | Editing sites called (of 30) | Germline (of 40) | Strand / junction (of 20 each) |
|---|---|---|---|---|---|
| Naive per-cell VAF ≥ 0.2, ≥ 2 alt reads | 0.17 | 0.70 | 30 | 40 | 20 / 20 |
| scMuta, no negative controls | 0.50 | 0.70 | 30 | 0.2 | 0 / 0 |
| **scMuta + negative controls** | **0.997** | 0.62 | 0.4 | 0 | 0 / 0 |

Artefact flags remove germline, strand and junction artefacts; **RNA editing is only removed when negative
controls exist**. Recall is limited by coverage and dropout (truth = mutant cells with ≥ 2 reads), and the
calibrated threshold trades some recall for precision (`--target-fp`).

## Usage

```bash
pip install -e .
# plate data (one BAM per cell)
scmuta call --plate cells.tsv --fasta GRCh38.fa --regions exons.bed --negative control_cells.txt --out out/
# droplet data
scmuta call --bam possorted.bam --tag CB --cells barcodes.tsv --stranded --fasta GRCh38.fa \
            --regions exons.bed --negative control_cells.txt --out out/
# known sites (e.g. an engineered mutation) with per-cell truth labels
scmuta call --plate cells.tsv --fasta GRCh38.fa --sites known.tsv --truth labels.csv --out out/
```

Outputs: `sites.tsv` (flags, error, cells called), `cell_calls.tsv` (cell, site, reads, posterior),
`calibration.tsv` (threshold vs held-out negative false-call rate), `summary.json`.

## Roadmap

- [x] First real-data truth test: FUS P525L knock-in Smart-seq2 (49/50 carriers, 47/48 controls; found the silent
      Y526 blocking edit blind) - `examples/fus_p525l.py`
- [ ] Stray-read (index hopping / cross-contamination) rate per site estimated from negatives
- [ ] More ground truth:
      public cell-line mixtures (HapMap / SMaHT), GoT / TARGET-seq data
- [ ] Germline heterozygous SNPs for allele-specific expression and dropout per cell (SCAN2-style)
- [ ] Editing modelled as a probability (cluster-wide, A>G context) rather than only screened
- [ ] Panel of normals across samples; Alu / repeat covariates
- [ ] Doublet and ambient-RNA terms; clone-level calling (phylogeny-aware)
- [ ] Speed: parallel per region; sparse matrices for large droplet data
- [ ] Benchmark vs SComatic, Monopogen, scVar, cellsnp-lite + VAF thresholds

## Licence

MIT
