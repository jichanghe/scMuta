# scMuta

**Per-cell somatic SNV calls from single-cell RNA-seq, with a false-call rate you can trust.**

[![tests](https://github.com/jichanghe/scMuta/actions/workflows/tests.yml/badge.svg)](https://github.com/jichanghe/scMuta/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23543571.svg)](https://doi.org/10.5281/zenodo.23543571)

Calling mutations in single cells from RNA is dominated by false signals: sequencing and reverse-transcription
errors, RNA editing, alignment artefacts at read ends and splice junctions, germline SNPs, and allelic dropout.
Most tools return hard calls without saying how many are wrong. scMuta returns a genotype for every cell at
every candidate site (**mutant / wild type / uncertain / not genotypable**) and measures its own false-call rate
on control cells that cannot carry the variants (another donor, another cell line, unedited cells). It works on
plate data (Smart-seq2, SMART-seq+, TARGET-seq; one BAM per cell) and droplet data (10x scRNA / snRNA; one BAM
with a cell-barcode tag).

Status: **0.2.0**, under active development. Interfaces may change before 1.0.

## How it works

| Step | What scMuta does |
|---|---|
| Candidates | Pooled pileup over cells; reads filtered by MAPQ, base quality, duplicates, mismatch load |
| Molecules | `--dedup umi` (droplet) or `start` (plate): counts molecules, not PCR copies |
| Site artefact flags | Alt reads compared with **reference reads at the same site**: strand bias, read-end and soft-clip enrichment (Fisher, one-sided); mismatches within 3 nt of a splice junction; homopolymer extension (run >= 4); base quality; known RNA-editing sites and editing clusters; germline-like sites |
| Per-cell genotype | Per-site mixture fitted by EM: wild type = Binomial(error), mutant = Beta-Binomial with a mutant allele fraction learned per site; error shrunk to a per-substitution background |
| Four states | mutant; wild type (the cell's reads favour wild type by >= 10:1); uncertain; not genotypable |
| **Calibration** | Leave-one-out on the control cells: each control cell's calls are counted only at sites the *other* controls do not carry and that would be reported. The posterior threshold is the lowest whose control false-call rate is <= `--target-fp`; the rate is reported |
| Reporting | A site is reported when called in >= 2 cells of interest; per-cell calls, the calibration table and a split-read reproducibility check are written out |

## Results on real data

**Engineered mutation, FUS P525L knock-in** (Smart-seq2, 50 P525L + 48 control cells): 49/50 carriers called
mutant, 48/48 controls called wild type. Without being told where to look, scMuta also found the silent T>C
blocking edit of the knock-in (chr16:31,191,435) in the same 49 cells. The one carrier called wild type has 0
mutant reads of 7 (likely a mislabelled cell).

**Driver mutation, JAK2 V617F** (TARGET-seq / SMART-seq+ cell lines, Rodriguez-Meira et al. 2019, GSE105451):
57 SET2 cells called mutant; 0 of 32 Jurkat and HSPC control cells.

**Genome-wide benchmark** (same dataset; `examples/cellline_benchmark.py`). Cells of each line and protocol are
split at random: one half defines the truth (variants private to SET2 or Jurkat - their germline SNPs and somatic
mutations - and sites no line carries), the other half is called. Case cells are SET2 + Jurkat mixed; HSPC cells
are the controls. 4,212 line-private variants, 76,942 covered carrier cell-site pairs, 733,907 sites that must stay
silent.

| Method | Sensitivity | Precision | False calls | False variant sites | Somatic calls at sites the controls carry |
|---|---|---|---|---|---|
| **scMuta** (threshold chosen by calibration, 0.9) | 59.0% | **97.3%** | 1,254 | 531 | 5,427 |
| scMuta at 0.999 | 57.5% | **98.7%** | - | **239** | - |
| >= 2 alt reads, VAF >= 0.1 | 63.9% | 96.8% | 1,651 | 718 | 114,730 |
| same, sites seen in controls removed | 57.3% | 96.8% | 1,465 | 638 | 8,758 |
| >= 1 alt read | 77.4% | 10.1% | 530,008 | 213,118 | 162,818 |

* **Calibration holds on independent truth**: false-positive rate estimated from the controls 5.4e-4, observed
  in the test cells 4.4e-4.
* At matched sensitivity scMuta makes ~2.5x fewer false calls and reports ~2.7x fewer false sites than the
  best threshold rule; results are the same on both protocols (SMART-seq+ / TARGET-seq: sensitivity 59.9% /
  58.4%, precision 97.0% / 97.5%).
* Sensitivity follows read depth (18% at 1-2 reads, 51% at 3-5, 68% at 6-10, 84% at > 10): with one or two
  reads a heterozygous mutation often shows no alt read. 7.3% of covered carrier pairs are called wild type -
  mostly cells expressing only the reference allele; "wild type" means *the reads support wild type*.

**Compared with other callers** on the same test cells and the same truth (`examples/competitors.py`).
SComatic was given the cell lines as cell types, which its model requires; the others were not. scAllele was run
on chromosomes 1 and 19 only, because a genome-wide run takes about an hour per cell, and all methods are scored
on those two chromosomes in that row.

| Method | Sensitivity | Precision | False calls | False variant sites | Runtime |
|---|---|---|---|---|---|
| **scMuta** | **58.3%** | **97.3%** | 1,254 | 531 | **3.2 h** |
| SComatic | 20.8% | 96.9% | **513** | **183** | 27.5 h |
| >= 2 alt reads, VAF >= 0.1 | 63.1% | 96.8% | 1,651 | 718 | - |
| same, sites seen in controls removed | 56.6% | 96.8% | 1,465 | 638 | - |

On chromosomes 1 and 19, where scAllele could also be run (15,897 carrier pairs):

| Method | Sensitivity | Precision | False calls | False variant sites |
|---|---|---|---|---|
| **scMuta** | **58.4%** | **96.8%** | 309 | 133 |
| SComatic | 21.0% | 96.5% | **121** | **39** |
| scAllele | 45.9% | 59.6% | 4,932 | 4,876 |
| >= 2 alt reads, VAF >= 0.1 | 62.2% | 96.4% | 366 | 159 |
| same, sites seen in controls removed | 55.7% | 96.4% | 332 | 144 |

The three fail in different directions. SComatic is conservative: high precision, but it misses four of five
carrier pairs, and its sensitivity *falls* at high depth (26% above 10 reads, against 84% for scMuta) because its
panel of normals and noisy-site filters remove the common variants that depth makes visible. scAllele recovers
variants wherever depth allows (81% above 10 reads) but calls germline and somatic alike, so two of every five
of its calls are wrong. Restricting the truth set to sites outside SComatic's panel of normals - the fair ground
for somatic callers - scMuta still made 1.7x more correct calls than SComatic (27,660 against 16,199) at
comparable precision (96.1% against 96.9%).

**Where RNA runs out: 3' snRNA of brain** (motor cortex, GSE174332: 3 C9orf72 FTLD donors and 3 controls; 42,478
nuclei, 22 autosomes and chrX, 3,031 Mb; `examples/snrna_pilot.py`, `snrna_summary.py`, `snrna_candidates.py`,
`snrna_site_check.py`, `snrna_classify.py`). 141 sites passed the somatic criteria, 7-61 per donor in proportion to
each donor's interrogable sequence (r = 0.88); none alters a protein. A read-level test of the donor's *non-carrier*
nuclei shows that **87 (62%) are inherited variants** (median 39% of their molecules carry the allele) that
nucleus-count germline rules miss in data this sparse; 28 are artefacts or weak, 17 unresolved, and **9 remain
candidates** (3 FTLD, 6 control; at most 2.6% of non-carrier molecules), the best an ARHGDIG 3'UTR variant in 8
excitatory neurons of one control. The calibration also states what could have been found. Median genotypable
positions per nucleus per chromosome (>= 2 molecules) run from 3,700 in excitatory neurons to 140 in microglia, a
**26-fold spread**:

| Cell type | Nuclei | Positions/nucleus | Mb detectable, clone at 10% | at 100% |
|---|---|---|---|---|
| Excitatory neurons | 10,994 | 3,700 | 2.01 | 18.68 |
| Interneurons | 983 | 2,160 | 0.12 | 1.69 |
| Astrocytes | 3,718 | 770 | 0.11 | 2.11 |
| OPC | 1,713 | 560 | 0.02 | 0.78 |
| Endothelial | 439 | 500 | 0.01 | 0.40 |
| Oligodendrocytes | 22,817 | 410 | 0.34 | 6.56 |
| Microglia | 1,065 | 140 | 0.002 | 0.12 |

Two things follow. The number of calls here **cannot** measure clonal mosaicism: only 0.62% of the sequence examined
is interrogable even for a clone in *every* nucleus of the best-covered cell type, so the expected number of
detections is far below one. And because coverage differs 26-fold between cell types, **somatic burden compared
across cell types in 3' snRNA data tracks coverage unless divided by genotypable positions** - excitatory neurons
will rank highest whatever the biology. `snrna_summary.py` reports the per-cell-type detection limit for exactly
this reason.

Runtimes are from a shared workstation and are indicative only.

## Install

```bash
pip install git+https://github.com/jichanghe/scMuta.git
# or, for development
git clone https://github.com/jichanghe/scMuta.git && cd scMuta && pip install -e ".[test]" && pytest -q
```
Requires Python >= 3.9, numpy, scipy, pandas, pysam. Linux or macOS.

## Quick start

```bash
bash examples/quickstart.sh       # builds a tiny synthetic data set and runs scmuta call (< 1 min)
```

Plate data (one BAM per cell):
```bash
scmuta call --plate cells.tsv --fasta GRCh38.fa --regions expressed_exons.bed \
            --negative control_cells.txt --out out/
```
`cells.tsv`: `cell<TAB>path/to/cell.bam` (coordinate-sorted, indexed). `control_cells.txt`: one cell per line.

Droplet data (10x, stranded):
```bash
scmuta call --bam possorted_genome_bam.bam --tag CB --cells barcodes.tsv --stranded \
            --fasta GRCh38.fa --regions expressed_exons.bed --negative control_cells.txt --out out/
```
Controls for droplet data are typically cells of another donor (pooled or genotype-demultiplexed samples).
Runs above `--block-sites` candidates (default 20,000) use the memory-light sparse pipeline automatically.

Known sites only (e.g. an engineered mutation), skipping discovery:
```bash
scmuta call --plate cells.tsv --fasta GRCh38.fa --sites known.tsv --truth labels.csv --out out/
```

### Outputs

| File | Content |
|---|---|
| `sites.tsv` | every candidate allele (0-based `pos`): filters, pooled VAF, error rate, mutant VAF, cells called, control calls |
| `cell_calls.tsv` | cell, site, alt reads, depth, posterior for every mutant call at a reported site |
| `calibration.tsv` | posterior threshold vs false-call rate in the control cells |
| `summary.json` | counts, chosen threshold, estimated false-positive rate, filter tallies |

### Gene and protein change

Add `--gtf gencode.gtf.gz` to `scmuta call` to write `variants_annotated.tsv` and `cell_calls_annotated.tsv`, or
annotate any table with `chrom, pos (0-based), ref, alt` afterwards:
```bash
scmuta annotate --table out/sites.tsv --gtf gencode.v50.annotation.gtf.gz --fasta GRCh38.fa --out sites_annotated.tsv
```
Columns added: `gene`, `gene_type`, `transcript_id` (MANE Select, else Ensembl canonical, else the longest coding
transcript), `region` (CDS, 5UTR, 3UTR, noncoding_exon, splice_site, intron, exon_other_isoform, intergenic; a site outside the chosen transcript but inside another isoform of the gene is intron or exon_other_isoform), `consequence`
(missense, synonymous, stop_gained, stop_lost, start_lost, splice_site, utr, noncoding, other_isoform, intronic, intergenic),
`codon_change`, `protein_change` (e.g. p.V617F) and `other_genes`. No external tools are needed. For databases
(gnomAD, ClinVar, COSMIC), pass `variants_annotated.tsv` to Ensembl VEP.

### Key options

`--target-fp` (default 0.001) false-call rate allowed in controls; `--min-cell-alt` (2) alt reads a cell needs to
be called mutant; `--min-alt-cells` (2) cells a site needs to be reported; `--editing-bed` known editing sites
(e.g. REDIportal); `--dedup none|umi|start`; `--stranded` for 10x.

## Reproducing the paper's analyses

| Analysis | Script |
|---|---|
| FUS P525L knock-in | `examples/fus_p525l.py` |
| Known-site truth test (JAK2 V617F) | `examples/known_site_test.py` |
| Genome-wide cell-line benchmark | `examples/cellline_benchmark.py`, diagnostics `examples/bench_diagnose.py` |
| SComatic / scAllele comparison | `examples/competitors.py` |
| Multi-donor 10x snRNA (motor cortex, C9orf72 FTLD and controls, GSE174332) | `examples/snrna_pilot.py` |
| Aggregate those runs and compute the clone detection limit | `examples/snrna_summary.py` |
| Collect and annotate the somatic-supported sites | `examples/snrna_candidates.py`, then `scmuta annotate` |
| Read-level check per site (molecule counts by nucleus group, optional pyGenomeTracks plots) | `examples/snrna_site_check.py` |
| Final classification (inherited / candidate / unresolved / artefact; Supplementary Table S1) | `examples/snrna_classify.py` |
| Memory probe for `count()` on one region | `examples/mem_probe.py` |
| Simulation | `examples/benchmark_sim.py` |

## Limitations

* RNA sees only expressed alleles: monoallelic expression and dropout limit sensitivity, and 3' droplet data
  cover few positions per cell - 3,700 positions per nucleus in excitatory neurons against 140 in microglia, so
  genome-wide discovery in 3' snRNA is power-limited rather than method-limited (see above).
* With controls from another individual, a mutation present in every case cell cannot be told apart from a
  germline SNP of that individual (flagged `germline_like`; cell-level calls are still reported).
* Calibration needs control cells that cannot carry the variants.

## Citation

Ji C. scMuta: calibrated per-cell somatic mutation calling from single-cell RNA sequencing. (in preparation).
Software: Ji C. scMuta 0.2.0. Zenodo. DOI of the 0.2.0 release (to be added once Zenodo archives it). See `CITATION.cff`.

## Licence

MIT
