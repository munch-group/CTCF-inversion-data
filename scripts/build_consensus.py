#!/usr/bin/env python
"""Cluster pooled peak summits into the consensus CTCF site table (README.md §4-5)."""
import argparse
import pandas as pd
import ctcf_lib as L

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--manifest", required=True)
ap.add_argument("--chunks", nargs="+", required=True, help="per-chunk .npz files")
ap.add_argument("--slop", type=int, default=L.SLOP)
ap.add_argument("--out", required=True, help="destination .tsv.gz")
a = ap.parse_args()

manifest = pd.read_csv(a.manifest, sep="\t")
P = L.load_chunks(a.chunks)
print(f"{P['summit'].size:,} peaks pooled from {len(a.chunks)} chunks")

sites = L.build_consensus(P, len(manifest), len(L.biosample_codes(manifest)), a.slop)
print(f"{len(sites):,} consensus CTCF sites")
print(f"  >=10 experiments : {(sites.n_experiments >= 10).sum():,}")
print(f"  >=100 experiments: {(sites.n_experiments >= 100).sum():,}")
print(f"  genome covered   : {sites.width.sum()/1e6:.1f} Mb")

sites.to_csv(a.out, sep="\t", index=False)
print(f"wrote {a.out}")
