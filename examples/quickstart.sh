#!/usr/bin/env bash
# Two-minute end-to-end example on a tiny synthetic data set (no downloads):
# 24 cells, one BAM each; 8 carry a somatic A>T, 8 are wild type, 8 are negative-control cells.
# The data also contain a splice-junction artefact and an RNA-editing-like site present in every cell.
set -euo pipefail
D=${1:-scmuta_quickstart}
python - "$D" <<'PY'
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath("examples/quickstart.sh")), "..", "tests"))
sys.path.insert(0, "tests")
from make_bams import build
build(sys.argv[1])
PY
scmuta call --plate "$D/cells.tsv" --fasta "$D/ref.fa" --regions "$D/regions.bed" \
            --negative "$D/negative.txt" --out "$D/out" > "$D/summary.txt"
echo "sites (with artefact flags):"; cut -f1-4,$(head -1 "$D/out/sites.tsv" | tr '\t' '\n' | grep -n -x filter | cut -d: -f1) "$D/out/sites.tsv"
echo; echo "cells called mutant:"; cut -f1-7 "$D/out/cell_calls.tsv"
echo; echo "expected: one variant (chrT:501 A>T) called in mut0..mut7 only; the junction artefact flagged"
echo "splice_junction and the editing-like site flagged in_negative_controls"
