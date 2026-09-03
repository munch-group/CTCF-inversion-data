#!/usr/bin/env python
"""Check the awk reimplementation reproduces the Python consensus (README.md §7)."""
import argparse
import sys

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--awk", required=True, help="output of merge_sites.awk")
ap.add_argument("--sites", required=True)
ap.add_argument("--chrom", default="chr21")
ap.add_argument("--out", required=True)
a = ap.parse_args()

cols = ["chrom", "start", "end", "summit", "n_experiments", "n_biosamples",
        "n_peaks", "max_signal", "mean_signal", "max_neglog10q"]
awk = pd.read_csv(a.awk, sep="\t", header=None, names=cols)
py = pd.read_csv(a.sites, sep="\t")
py = py[py.chrom == a.chrom].reset_index(drop=True)

lines = [f"awk reimplementation vs python consensus on {a.chrom}",
         f"  awk sites    {len(awk):,}",
         f"  python sites {len(py):,}"]
ok = len(awk) == len(py)
if ok:
    for c in ["start", "end", "summit", "n_experiments", "n_biosamples", "n_peaks"]:
        bad = int((awk[c].to_numpy() != py[c].to_numpy()).sum())
        ok &= bad == 0
        lines.append(f"  {c:<16} exact mismatches: {bad}")
    for c in ["max_signal", "mean_signal", "max_neglog10q"]:
        d = float(np.abs(awk[c].to_numpy() - py[c].to_numpy()).max())
        ok &= d < 1e-3
        lines.append(f"  {c:<16} max abs diff: {d:.3e}")
else:
    lines.append("  ROW COUNT MISMATCH")
lines.append(f"\nRECIPE REPRODUCES PYTHON OUTPUT: {ok}")

report = "\n".join(lines)
print(report)
with open(a.out, "w") as fh:
    fh.write(report + "\n")
sys.exit(0 if ok else 1)
