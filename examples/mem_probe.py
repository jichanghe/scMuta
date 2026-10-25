#!/usr/bin/env python
"""Memory probe for count(): one donor, one region, two versions of scmuta/pileup.py, each in its own process.

  python mem_probe.py --bam STARsolo/210/Aligned.sortedByCoord.out.bam --donor 210 \
      --meta ~/snrna_genome/chr12/coverage_per_nucleus.tsv --fasta GRCh38.fa --region chr11:60000000-70000000 \
      --pkg old=/path/to/old_scMuta --pkg new=/path/to/new_scMuta

Discovers candidates in the region once (with the first package), then times count() for each package and prints
peak resident memory, run time, number of fetch windows and totals (which must agree between versions).
"""
import argparse
import json
import os
import subprocess
import sys

CHILD = r'''
import sys, json, time, resource, pandas as pd
sys.path.insert(0, sys.argv[1])
from scmuta.pileup import count, ReadFilter
import scmuta.pileup as P
args = json.loads(sys.argv[2])
sites = pd.read_csv(args["sites"], sep="\t")
keep = set(open(args["keep"]).read().split())
t = time.time()
kw = dict(umi_tag="UB", dedup="umi")
import inspect
if "sparse_out" in inspect.signature(count).parameters:
    kw["sparse_out"] = True
A, D, s, cells = count(args["bam"], args["fasta"], sites, ReadFilter(), "CB", keep, **kw)
peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024
mw = inspect.signature(count).parameters.get("max_window")
print(json.dumps({"peak_MB": peak, "seconds": round(time.time() - t, 1), "cells": len(cells),
                  "alt_total": int(A.sum()), "depth_total": int(D.sum()),
                  "max_window": None if mw is None else mw.default, "file": P.__file__}))
'''


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--donor", required=True)
    ap.add_argument("--meta", required=True, help="coverage_per_nucleus.tsv from any chromosome run (barcodes per donor)")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--region", default="chr11:60000000-70000000")
    ap.add_argument("--pkg", action="append", required=True, help="NAME=folder containing the scmuta package")
    ap.add_argument("--work", default="mem_probe_work")
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    import pandas as pd
    m = pd.read_csv(os.path.expanduser(a.meta), sep="\t", dtype={"donor": str})
    keep = m[m.donor == a.donor].barcode.tolist()
    keep_f = os.path.join(a.work, "keep.txt")
    open(keep_f, "w").write("\n".join(keep))
    chrom, rng = a.region.split(":")
    s0, s1 = (int(x.replace(",", "")) for x in rng.split("-"))
    sites_f = os.path.join(a.work, "sites.tsv")
    first = a.pkg[0].split("=", 1)[1]
    if not os.path.exists(sites_f):
        sys.path.insert(0, first)
        from scmuta.pileup import discover, ReadFilter
        print(f"discovering candidates in {a.region} for donor {a.donor} ({len(keep)} nuclei) ...", flush=True)
        regions = [(chrom, x, min(s1, x + 5_000_000)) for x in range(s0, s1, 5_000_000)]
        cd = discover(a.bam, a.fasta, regions, ReadFilter(), "CB", set(keep), min_alt_reads=3, min_alt_cells=2)
        cd[["chrom", "pos", "ref", "alt"]].to_csv(sites_f, sep="\t", index=False)
    sites = pd.read_csv(sites_f, sep="\t")
    p = sites.pos.sort_values().to_numpy()
    gaps = int((p[1:] - p[:-1] > 5000).sum())
    print(f"{len(sites)} candidate sites; spacing median {int(pd.Series(p).diff().median())} nt; "
          f"windows: {gaps + 1} with the 5 kb rule alone", flush=True)
    args = json.dumps({"sites": sites_f, "keep": keep_f, "bam": a.bam, "fasta": a.fasta})
    for spec in a.pkg:
        name, folder = spec.split("=", 1)
        print(f"\n[{name}] counting ...", flush=True)
        r = subprocess.run([sys.executable, "-c", CHILD, folder, args], capture_output=True, text=True)
        out = [l for l in r.stdout.splitlines() if l.startswith("{")]
        if r.returncode != 0 or not out:
            print(f"[{name}] FAILED (exit {r.returncode})\n{r.stderr[-2000:]}")
            continue
        print(f"[{name}] " + out[-1], flush=True)


if __name__ == "__main__":
    main()
