#!/usr/bin/env python
"""Download and parse one chunk of the ENCODE peak files (README.md §3).

The manifest is split by row index modulo --chunks so the chunk count is fixed
at workflow-definition time; experiment and biosample indices are derived from
the manifest, so they are consistent across chunks.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import ctcf_lib as L

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--manifest", required=True)
ap.add_argument("--chunk", type=int, required=True)
ap.add_argument("--chunks", type=int, required=True)
ap.add_argument("--workers", type=int, default=4)
ap.add_argument("--out", required=True, help="destination .npz")
a = ap.parse_args()

manifest = pd.read_csv(a.manifest, sep="\t")
bio_codes = L.biosample_codes(manifest)
mine = manifest.iloc[a.chunk::a.chunks]
print(f"chunk {a.chunk}/{a.chunks}: {len(mine)} experiments")


def work(row):
    code, summit, signal, negq = L.parse_peaks(L.fetch(row.url))
    n = summit.size
    return (code, summit, signal, negq,
            np.full(n, row.Index, dtype=np.int16),
            np.full(n, bio_codes[row.biosample], dtype=np.int16))


with ThreadPoolExecutor(a.workers) as pool:
    parts = list(pool.map(work, mine.itertuples()))

keys = ("chrom", "summit", "signal", "negq", "exp", "bio")
arrays = {k: np.concatenate([p[i] for p in parts]) for i, k in enumerate(keys)}
print(f"  {arrays['summit'].size:,} peaks")
np.savez_compressed(a.out, **arrays)
print(f"wrote {a.out}")
