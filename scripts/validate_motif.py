#!/usr/bin/env python
"""Independent check that the sites are real CTCF sites in hg38 (METHODS.md §9).

Scans the JASPAR CTCF motif MA0139.1 around every site summit on chr21 and
compares the hit rate against random chr21 positions and against the same
summits shifted 5 kb. Sequence evidence never entered the compilation, so a
recurrence-dependent enrichment that collapses under the shift confirms both
the assembly and the biological identity of the sites. Wrong-assembly
coordinates would sit at the background rate throughout.
"""
import argparse
import gzip
import json
import sys

import numpy as np
import pandas as pd
import ctcf_lib as L

REL_CUTOFF = 0.80
FLANK = 100

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--sites", required=True)
ap.add_argument("--chr21", required=True, help="hg38 chr21.fa.gz")
ap.add_argument("--pfm", required=True, help="JASPAR MA0139.1 JSON")
ap.add_argument("--out", required=True)
a = ap.parse_args()

with gzip.open(a.chr21, "rt") as fh:
    seq = "".join(l.strip() for l in fh if not l.startswith(">"))
S = L.encode_sequence(seq)
lo, rc, mn, mx = L.load_pwm(json.load(open(a.pfm))["pfm"])
score = lambda c: L.best_relscore(S, c, lo, rc, mn, mx, FLANK)

sites = pd.read_csv(a.sites, sep="\t")
c21 = sites[sites.chrom == "chr21"]
c21 = c21[(c21.summit > 200) & (c21.summit < len(S) - 200)].copy()
c21["rel"] = score(c21.summit.to_numpy())

rng = np.random.default_rng(0)
ok_pos = np.flatnonzero(S >= 0)
ok_pos = ok_pos[(ok_pos > 200) & (ok_pos < len(S) - 200)]

rate = lambda x: (np.asarray(x) >= REL_CUTOFF).mean() * 100
fmt = lambda x: f"{rate(x):5.1f}%  (n={len(x):,})"

n = c21.n_experiments
bg = score(rng.choice(ok_pos, 5000, replace=False))
shifted = c21.loc[n >= 100, "summit"].to_numpy() + 5000
shifted = shifted[shifted < len(S) - 200]
sh = score(shifted)

lines = [f"CTCF motif within +/-{FLANK}bp of summit (MA0139.1, rel >= {REL_CUTOFF})",
         f"  random chr21 positions    : {fmt(bg)}",
         f"  sites in    1 experiment  : {fmt(c21.loc[n == 1, 'rel'])}",
         f"  sites in  2-9 experiments : {fmt(c21.loc[n.between(2, 9), 'rel'])}",
         f"  sites in >=10 experiments : {fmt(c21.loc[n >= 10, 'rel'])}",
         f"  sites in >=100 experiments: {fmt(c21.loc[n >= 100, 'rel'])}",
         f"  sites in >=300 experiments: {fmt(c21.loc[n >= 300, 'rel'])}",
         f"  ... same sites, +5kb      : {fmt(sh)} <- negative control"]

checks = {
    "recurrent sites (>=100 expts) are >80% motif-positive": rate(c21.loc[n >= 100, "rel"]) > 80,
    "random background is <15%": rate(bg) < 15,
    "+5kb shifted control is <20%": rate(sh) < 20,
    "motif rate rises with recurrence": (rate(c21.loc[n == 1, "rel"])
                                         < rate(c21.loc[n.between(2, 9), "rel"])
                                         < rate(c21.loc[n >= 10, "rel"])
                                         < rate(c21.loc[n >= 100, "rel"])),
}
lines.append("")
for k, v in checks.items():
    lines.append(f"  {'PASS' if v else 'FAIL'}  {k}")
lines.append(f"\nMOTIF VALIDATION: {'PASS' if all(checks.values()) else 'FAIL'}")

report = "\n".join(lines)
print(report)
with open(a.out, "w") as fh:
    fh.write(report + "\n")
sys.exit(0 if all(checks.values()) else 1)
