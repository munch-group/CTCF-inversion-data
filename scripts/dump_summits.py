#!/usr/bin/env python
"""Write the pooled peak summits for one chromosome as an annotated BED.

Feeds the dependency-free awk reimplementation of the clustering step, so the
consensus algorithm can be checked against a second implementation.
"""
import argparse

import numpy as np
import pandas as pd
import ctcf_lib as L

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--manifest", required=True)
ap.add_argument("--chunks", nargs="+", required=True)
ap.add_argument("--chrom", default="chr21")
ap.add_argument("--out", required=True)
a = ap.parse_args()

manifest = pd.read_csv(a.manifest, sep="\t")
code = L.CHROM_CODE[a.chrom]
P = L.load_chunks(a.chunks)

m = P["chrom"] == code
exp_name = manifest["experiment"].to_numpy()
bio_name = manifest["biosample"].to_numpy()
inv_bio = {v: k for k, v in L.biosample_codes(manifest).items()}

with open(a.out, "w") as fh:
    for s, sig, q, e, b in zip(P["summit"][m], P["signal"][m], P["negq"][m],
                               P["exp"][m], P["bio"][m]):
        fh.write(f"{a.chrom}\t{s}\t{s+1}\t{exp_name[e]}\t{inv_bio[b]}\t{sig:.6f}\t{q:.6f}\n")
print(f"wrote {int(m.sum()):,} summits to {a.out}")
