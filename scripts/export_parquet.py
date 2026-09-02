#!/usr/bin/env python
"""Convert a site table into a parquet dataset of GitHub-sized part files.

`pd_lfs.write_parquet` shards the frame into `part-*.parquet` files each below
`--max-bytes` and writes a `_manifest.json` index alongside them, which lets
`pd_lfs.read_parquet` read the dataset over plain HTTPS without directory
listing. Column dtypes are downcast on the way out and recorded in the
manifest, so reading restores the original dtypes transparently.

Conversion is lossless: floats are left at full precision (`float_decimals` is
deliberately not set), so the round-trip reproduces the source frame exactly.
"""
import argparse
import os

from pd_lfs import write_parquet

import ctcf_lib as L

# GitHub warns above 50 MB and blocks above 100 MB; leave headroom under the warning.
MAX_BYTES = 40_000_000

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--sites", required=True, help="source .tsv.gz")
ap.add_argument("--out", required=True, help="destination dataset directory")
ap.add_argument("--max-bytes", type=int, default=MAX_BYTES)
a = ap.parse_args()

# load through the canonical loader so the dataset carries the compact
# dtypes (int32/int16/float32/category) rather than pandas' CSV defaults
df = L.load_ctcf_sites(a.sites)
write_parquet(df, a.out, max_bytes=a.max_bytes)

parts = sorted(f for f in os.listdir(a.out) if f.endswith(".parquet"))
sizes = [os.path.getsize(os.path.join(a.out, f)) for f in parts]
print(f"{a.sites} -> {a.out}")
print(f"  {len(df):,} rows x {df.shape[1]} columns")
print(f"  {len(parts)} part file(s), {sum(sizes)/1e6:.1f} MB total")
for f, s in zip(parts, sizes):
    print(f"    {s/1e6:8.2f} MB  {f}")

over = [f for f, s in zip(parts, sizes) if s > a.max_bytes]
assert not over, f"part files exceed --max-bytes: {over}"
print(f"  all part files below {a.max_bytes/1e6:.0f} MB")
