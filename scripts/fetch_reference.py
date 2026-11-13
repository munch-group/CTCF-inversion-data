#!/usr/bin/env python
"""Download one external reference file used by the motif, liftover or validation steps."""
import argparse
import ctcf_lib as L

URLS = {
    "chain": L.CHAIN_URL,                # UCSC hg38 -> hs1 (CHM13v2.0) liftOver chains
    "jaspar": L.JASPAR_URL,              # JASPAR CTCF matrix MA0139.1
    "jaspar_track": L.JASPAR_TRACK_URL,  # JASPAR's genome-wide MA0139.1 hg38 predictions
    "hs1_sizes": L.HS1_SIZES_URL,        # CHM13 chromosome lengths, for bounds checking
    # hg38 chromosome sequences, for the motif scan (chr21 also for the validation)
    **{c: L.HG38_CHROM_URL.format(chrom=c) for c in L.CHROMS},
}

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--what", required=True, choices=sorted(URLS))
ap.add_argument("--out", required=True)
a = ap.parse_args()

with open(a.out, "wb") as fh:
    fh.write(L.fetch(URLS[a.what]))
print(f"wrote {a.out} from {URLS[a.what]}")
