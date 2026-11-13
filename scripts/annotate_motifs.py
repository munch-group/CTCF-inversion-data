#!/usr/bin/env python
"""Merge the per-chromosome motif calls into the consensus site table (README.md §9.2)."""
import argparse

import pandas as pd
import ctcf_lib as L

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--sites", required=True, help="consensus site table (.tsv.gz)")
ap.add_argument("--scans", nargs="+", required=True, help="per-chromosome scan_motifs.py outputs")
ap.add_argument("--out", required=True, help="destination .tsv.gz")
a = ap.parse_args()

sites = pd.read_csv(a.sites, sep="\t")
scans = pd.concat([pd.read_csv(p, sep="\t", compression="gzip") for p in a.scans],
                  ignore_index=True)
assert len(scans) == len(sites), f"{len(scans):,} scanned rows for {len(sites):,} sites"

# (chrom, start) is unique: sites on a chromosome are disjoint intervals
out = sites.merge(scans, on=["chrom", "start"], how="left", validate="one_to_one")
assert len(out) == len(sites) and out.motif_strand.notna().all(), "sites without a motif call"
out = out.sort_values(["chrom", "start"], key=lambda s: s.map(L.CHROM_CODE) if s.name == "chrom" else s,
                      ignore_index=True)

n = len(out)
has = out.motif_strand != "."
plus = (out.motif_strand == "+").sum()
print(f"{n:,} sites, {has.sum():,} oriented ({has.mean()*100:.1f}%), "
      f"{plus / max(has.sum(), 1) * 100:.1f}% of those on the plus strand")
for lo, hi, lab in [(1, 1, "1"), (2, 9, "2-9"), (10, 99, "10-99"), (100, 10**6, ">=100")]:
    m = out.n_experiments.between(lo, hi)
    print(f"  {lab:>6} experiments: {has[m].mean()*100:5.1f}% oriented  (n={m.sum():,})")
both = (out.n_motifs_plus > 0) & (out.n_motifs_minus > 0)
print(f"  motifs in both orientations within the site: {both.sum():,} ({both.mean()*100:.2f}%)")

out.to_csv(a.out, sep="\t", index=False)
print(f"wrote {a.out}")
