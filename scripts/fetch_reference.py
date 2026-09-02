#!/usr/bin/env python
"""Download one external reference file used by the liftover or validation steps."""
import argparse
import ctcf_lib as L

URLS = {
    "chain": L.CHAIN_URL,          # UCSC hg38 -> hs1 (CHM13v2.0) liftOver chains
    "chr21": L.CHR21_URL,          # hg38 chr21 sequence, for the motif validation
    "jaspar": L.JASPAR_URL,        # JASPAR CTCF matrix MA0139.1
    "hs1_sizes": L.HS1_SIZES_URL,  # CHM13 chromosome lengths, for bounds checking
}

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--what", required=True, choices=sorted(URLS))
ap.add_argument("--out", required=True)
a = ap.parse_args()

with open(a.out, "wb") as fh:
    fh.write(L.fetch(URLS[a.what]))
print(f"wrote {a.out} from {URLS[a.what]}")
