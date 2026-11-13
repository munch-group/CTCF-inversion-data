#!/usr/bin/env python
"""Check the motif orientation calls against JASPAR's own hg38 predictions (README.md §9.3).

JASPAR publishes genome-wide MA0139.1 hits for hg38 (relative score >= 0.8,
scanned with PWMScan). Those predictions never entered the compilation, so
agreement in position and strand between our calls and theirs confirms both the
scan and the strand convention. Small differences in PWM construction make a
few hits near the 0.8 cutoff appear in one set only, so the checks are rates,
not identities.
"""
import argparse
import sys

import numpy as np
import pandas as pd
import ctcf_lib as L

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--sites", required=True, help="annotated site table (.tsv.gz)")
ap.add_argument("--track", required=True, help="JASPAR MA0139.1 hg38 .tsv.gz")
ap.add_argument("--out", required=True)
a = ap.parse_args()

sites = pd.read_csv(a.sites, sep="\t")
track = pd.read_csv(a.track, sep="\t", header=None, names=L.JASPAR_TRACK_COLUMNS)
track = track[track.chrom.isin(L.CHROMS)]
W = int((track.end - track.start).iloc[0])

called = sites[sites.motif_strand != "."]
key = lambda df, start, strand: (df.chrom.astype(str) + ":" + df[start].astype(str) + df[strand].astype(str))
jas = set(key(track, "start", "strand"))
exact = key(called, "motif_start", "motif_strand").isin(jas).to_numpy()

# same position but opposite strand would mean a strand-convention error
flipped = called.motif_strand.map({"+": "-", "-": "+"})
opp = (called.chrom.astype(str) + ":" + called.motif_start.astype(str) + flipped).isin(jas).to_numpy()

# JASPAR's PWM is built with different pseudocounts, so its relative scores run
# ~0.01 lower than ours and hits just above 0.80 drop out of its track; above
# a small margin the two scans must agree
rel = called.motif_rel.to_numpy()
EDGE = L.MOTIF_REL_CUTOFF + 0.02
clear = rel >= EDGE

# the converse: a JASPAR hit within the summit window of a recurrent site must
# lie in a site that got oriented (it is either our call or a weaker second
# motif in the same window)
core = sites[sites.n_experiments >= 100]
in_window = np.zeros(len(track), bool)
oriented = np.zeros(len(track), bool)
for c, g in track.groupby("chrom", observed=True):
    cs = core[core.chrom == c]
    if cs.empty:
        continue
    order = np.argsort(cs.summit.to_numpy())
    s = cs.summit.to_numpy(np.int64)[order]
    has = (cs.motif_strand != ".").to_numpy()[order]
    p = g.start.to_numpy(np.int64)
    j = np.searchsorted(s, p)
    for k in (np.maximum(j - 1, 0), np.minimum(j, s.size - 1)):
        near = (p >= s[k] - L.MOTIF_FLANK) & (p <= s[k] + L.MOTIF_FLANK - W + 1)
        rows = track.index.get_indexer(g.index)
        in_window[rows] |= near
        oriented[rows] |= near & has[k]
core_hits_oriented = oriented[in_window].mean() * 100

plus_frac = (called.motif_strand == "+").mean() * 100
n_core = (core.motif_strand != ".").mean() * 100

lines = [f"motif orientation vs JASPAR 2022 hg38 track for MA0139.1 ({len(track):,} hits on chr1-22,X,Y)",
         f"  sites with an oriented motif                 : {len(called):,} / {len(sites):,}",
         f"  same position and strand in JASPAR           : {exact.mean()*100:5.1f}%",
         f"    ... among calls with motif_rel >= {EDGE:.2f}    : {exact[clear].mean()*100:5.1f}%  (n={clear.sum():,})",
         f"    ... among calls with motif_rel <  {EDGE:.2f}    : {exact[~clear].mean()*100:5.1f}%  (n={(~clear).sum():,})",
         f"  same position, opposite strand               : {opp.mean()*100:5.1f}%  <- must be ~0",
         f"  oriented sites on the plus strand            : {plus_frac:5.1f}%",
         f"  >=100-experiment sites oriented              : {n_core:5.1f}%",
         f"  JASPAR hits at >=100-expt summits, site oriented: {core_hits_oriented:5.1f}%  (n={in_window.sum():,})"]

checks = {
    f"calls with motif_rel >= {EDGE:.2f} agree with JASPAR in position and strand (>99%)": exact[clear].mean() > 0.99,
    "no called motif has the opposite strand in JASPAR (<0.1%)": opp.mean() < 0.001,
    "strand balance is 45-55% plus": 45 < plus_frac < 55,
    ">80% of >=100-experiment sites are oriented": n_core > 80,
    "every JASPAR hit at a recurrent summit lies in an oriented site (>99%)": core_hits_oriented > 99,
}
lines.append("")
for k, v in checks.items():
    lines.append(f"  {'PASS' if v else 'FAIL'}  {k}")
lines.append(f"\nORIENTATION VALIDATION: {'PASS' if all(checks.values()) else 'FAIL'}")

report = "\n".join(lines)
print(report)
with open(a.out, "w") as fh:
    fh.write(report + "\n")
sys.exit(0 if all(checks.values()) else 1)
