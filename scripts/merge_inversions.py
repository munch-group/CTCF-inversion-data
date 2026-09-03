#!/usr/bin/env python
"""Combine the per-species inversion tables into one catalogue."""
import argparse

import pandas as pd

import ape_lib as A

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--tables", nargs="+", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

df = pd.concat([pd.read_csv(p, sep="\t") for p in a.tables], ignore_index=True)
order = {s: i for i, s in enumerate(A.OUTGROUP_ORDER + ["human_CHM13"])}
df["_o"] = df.species.map(order).fillna(99)
df = df.sort_values(["_o", "ref_chrom", "ref_start"]).drop(columns="_o")

print("inversions per species (count, Mb of reference):")
g = df.groupby("species").agg(n=("ref_span", "size"), Mb=("ref_span", lambda x: x.sum() / 1e6))
print(g.round(1).to_string())
print(f"\ntotal rows: {len(df):,}")

df.to_csv(a.out, sep="\t", index=False)
print(f"wrote {a.out}")
