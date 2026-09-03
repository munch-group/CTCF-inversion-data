#!/usr/bin/env python
"""Consistency checks on the derived CHM13 coordinates (README.md §11.5-11.6).

The decisive check is monotonicity within a chain: sites ordered by hg38
position must lift in ascending order on a plus chain and descending on a minus
chain. Order is deliberately NOT checked across a whole chromosome, because
sites switch chains wherever the two references disagree on orientation -- that
is signal, not error, and the flipped blocks are reported below.
"""
import argparse
import sys

import numpy as np
import pandas as pd
import ctcf_lib as L

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--sites", required=True, help="dual-coordinate table")
ap.add_argument("--chain", required=True)
ap.add_argument("--chrom-sizes", required=True, help="hs1.chrom.sizes")
ap.add_argument("--out", required=True)
a = ap.parse_args()

sz = pd.read_csv(a.chrom_sizes, sep="\t", header=None, names=["c", "n"])
sizes = dict(zip(sz.c, sz.n))
d = pd.read_csv(a.sites, sep="\t")
d["chm13_chrom"] = d["chm13_chrom"].fillna("")
u = d[d.chm13_status == "unique"]

checks = {}
lifted = set(d.loc[d.chm13_chrom != "", "chm13_chrom"])
checks["all lifted chromosomes exist in hs1"] = lifted <= set(sizes)
lim = u.chm13_chrom.map(sizes)
checks["within CHM13 chromosome bounds"] = bool(((u.chm13_start >= 0) & (u.chm13_end <= lim)).all())
checks["summit inside lifted interval"] = bool(
    ((u.chm13_summit >= u.chm13_start) & (u.chm13_summit < u.chm13_end)).all())
checks["start < end"] = bool((u.chm13_start < u.chm13_end).all())
checks["unmapped rows blanked"] = bool(
    (d.loc[d.chm13_status == "unmapped", ["chm13_start", "chm13_end", "chm13_summit"]] == -1).all().all())
w = u.chm13_end - u.chm13_start
checks["width within tolerance for every unique site"] = bool(
    (np.abs(w - u.width) <= L.GAP_TOLERANCE).all())

# the decisive check
chains = L.load_chains(a.chain)
_, _, _, chain_id = L.lift(chains, d["chrom"].to_numpy(dtype=object),
                           d["summit"].to_numpy(np.int64))
d["_chain"] = chain_id
bad = checked = 0
for _, g in d[d.chm13_status == "unique"].groupby("_chain"):
    if len(g) < 2:
        continue
    checked += 1
    g = g.sort_values("summit")
    mono = (g.chm13_summit.is_monotonic_increasing if g.chm13_strand.iloc[0] == "+"
            else g.chm13_summit.is_monotonic_decreasing)
    bad += not mono
checks[f"within-chain monotonicity ({checked} chains, {bad} violating)"] = bad == 0

lines = ["liftover consistency checks", ""]
lines += [f"  {'PASS' if v else 'FAIL'}  {k}" for k, v in checks.items()]
lines += ["", f"  width identical for {(w == u.width).mean()*100:.2f}% of unique sites",
          "", "mapping status:"]
lines += [f"  {k:<9}{v:>8,}" for k, v in d.chm13_status.value_counts().items()]
lines += ["", f"  minus strand (reference orientation flip): {(d.chm13_strand=='-').sum():,}",
          f"  lifted to a different chromosome:          "
          f"{((d.chm13_chrom != '') & (d.chm13_chrom != d.chrom)).sum():,}"]

lines += ["", "contiguous reference-orientation-flipped blocks (>=20 sites):",
          "  these are loci where GRCh38 and CHM13 carry opposite alleles;",
          "  polarise every inversion callset to one reference before merging."]
m = d[d.chm13_strand == "-"]
rows = []
for c, g in m.groupby("chrom"):
    g = g.sort_values("start")
    for _, blk in g.groupby((g.start.diff() > 2_000_000).cumsum()):
        if len(blk) >= 20:
            rows.append((c, int(blk.start.min()), int(blk.end.max()),
                         (blk.end.max() - blk.start.min()) / 1e6, len(blk),
                         blk.chm13_chrom.iloc[0]))
for r in sorted(rows, key=lambda x: -x[3]):
    lines.append(f"  {r[0]:<6}{r[1]:>12,}-{r[2]:<12,} {r[3]:6.2f} Mb  {r[4]:>4} sites -> {r[5]}")

lines.append(f"\nLIFTOVER VALIDATION: {'PASS' if all(checks.values()) else 'FAIL'}")
report = "\n".join(lines)
print(report)
with open(a.out, "w") as fh:
    fh.write(report + "\n")
sys.exit(0 if all(checks.values()) else 1)
