#!/usr/bin/env python
"""Stream one all-to-all chain file and keep only its chain header lines.

The chain files are 13-32 MB gzipped each; the alignment blocks are not needed
to locate inversions, only the per-chain header (target/query span, orientation
and score). Streaming and discarding the blocks keeps ~1.4 GB of alignment off
disk while retaining everything the caller uses.
"""
import argparse
import sys
import zlib

import requests

import ape_lib as A

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--ref", required=True, help="target assembly, e.g. hg38")
ap.add_argument("--query", required=True, help="query assembly accession or hs1")
ap.add_argument("--out", required=True)
ap.add_argument("--timeout", type=int, default=900)
a = ap.parse_args()

url = A.chain_url(a.ref, a.query)
r = requests.get(url, stream=True, timeout=a.timeout)
r.raise_for_status()

dec = zlib.decompressobj(zlib.MAX_WBITS | 16)
tail = ""
n = 0
with open(a.out, "w") as fh:
    for chunk in r.iter_content(1 << 20):
        text = tail + dec.decompress(chunk).decode("ascii", "replace")
        lines = text.split("\n")
        tail = lines.pop()
        for line in lines:
            if line.startswith("chain"):
                fh.write(line + "\n")
                n += 1
    text = tail + dec.flush().decode("ascii", "replace")
    for line in text.split("\n"):
        if line.startswith("chain"):
            fh.write(line + "\n")
            n += 1

# a truncated transfer leaves the gzip stream unfinished; fail rather than
# silently writing a partial header set
if not dec.eof:
    sys.exit(f"ERROR: gzip stream from {url} ended prematurely")
print(f"{url}\n  {n:,} chain headers -> {a.out}")
