#!/usr/bin/env python
"""Download the ENCODE batch-metadata table for CTCF narrowPeak files on GRCh38.

One row per file; the columns used downstream are documented in METHODS.md §2.1.
"""
import argparse
import ctcf_lib as L

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--out", required=True, help="destination TSV")
a = ap.parse_args()

with open(a.out, "wb") as fh:
    fh.write(L.fetch(L.METADATA_URL))
print(f"wrote {a.out}")
