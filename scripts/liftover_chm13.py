#!/usr/bin/env python
"""Append T2T-CHM13v2.0 coordinates to the hg38 site table (README.md §11).

ENCODE has no CTCF ChIP-seq on CHM13, so CHM13 coordinates are derived by
lifting the native GRCh38 coordinates through the UCSC hg38->hs1 chains. Unlike
UCSC `liftOver` this keeps the alignment strand: a site on a minus-strand chain
sits in a segment whose orientation differs between the two references.
"""
import argparse

import numpy as np
import pandas as pd
import ctcf_lib as L

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--sites", required=True)
ap.add_argument("--chain", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

chains = L.load_chains(a.chain)
print(f"{len(chains):,} chains")

sites = pd.read_csv(a.sites, sep="\t")
chrom = sites["chrom"].to_numpy(dtype=object)

a_chrom, a_pos, _, a_chain = L.lift(chains, chrom, sites["start"].to_numpy(np.int64))
s_chrom, s_pos, s_strand, s_chain = L.lift(chains, chrom, sites["summit"].to_numpy(np.int64))
b_chrom, b_pos, _, b_chain = L.lift(chains, chrom, sites["end"].to_numpy(np.int64) - 1)

mapped = s_pos >= 0
same = mapped & (a_chain == s_chain) & (b_chain == s_chain) & (a_pos >= 0) & (b_pos >= 0)

lo = np.minimum(a_pos, b_pos)
hi = np.maximum(a_pos, b_pos)
out_start = np.where(same, lo, -1)
out_end = np.where(same, hi + 1, -1)

# Both ends in one chain but on opposite sides of a large chain gap gives a
# wildly wrong width; small indels are legitimate.
delta = np.abs((out_end - out_start) - sites["width"].to_numpy())
status = np.where(~mapped, "unmapped",
          np.where(same & (delta > L.GAP_TOLERANCE), "gapped",
           np.where(same, "unique", "partial")))

sites["chm13_chrom"] = np.where(mapped, s_chrom, "")
sites["chm13_start"] = out_start
sites["chm13_end"] = out_end
sites["chm13_summit"] = np.where(mapped, s_pos, -1)
sites["chm13_strand"] = np.where(mapped, np.where(s_strand < 0, "-", "+"), "")
sites["chm13_status"] = status

n = len(sites)
print(f"mapped   {mapped.sum():>7,} / {n:,} ({mapped.mean()*100:.2f}%)")
for st in ("unique", "gapped", "partial", "unmapped"):
    print(f"  {st:<9}{(status == st).sum():>7,}")
print(f"minus strand (reference orientation flip): {(sites.chm13_strand == '-').sum():,}")
print(f"lifted to a different chromosome:          {(mapped & (s_chrom != chrom)).sum():,}")

sites.to_csv(a.out, sep="\t", index=False)
print(f"wrote {a.out}")
