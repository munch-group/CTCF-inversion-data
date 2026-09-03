#!/usr/bin/env python
"""Call inversions between two assemblies from all-to-all chain headers.

A chain records the orientation in which a target segment aligns to the query,
so an inversion is a run of chains opposite to the surrounding synteny. Three
things make a naive "minus strand = inversion" wrong, and all are handled here:

1. **Assembly orientation conventions.** Whole ape chromosomes are frequently
   stored reverse-complemented relative to their human orthologue -- chimpanzee
   chr1 is one -- so the minus strand is the *normal* state there. The dominant
   orientation is established per orthologous chromosome pair (score-weighted)
   and inversions are called against it.

2. **Paralogy.** These are raw Cactus all-to-all chains, not netted chains, so
   every region also aligns to its paralogues; target spans sum to well over a
   genome length. Restricting to the single orthologous query chromosome per
   reference chromosome keeps the syntenic backbone.

3. **Spurious small chains.** Below ~1 kb, minority-orientation chains are
   dominated by alignment noise: two human assemblies carry 7,440 of them,
   against chimpanzee's 13,915, despite ~6 Myr less divergence. The colinearity
   test below is what separates signal from noise at small sizes.

**Two categories, two criteria.** Large pericentric inversions and small local
inversions need different tests, and a single compromise criterion serves both
badly:

*small* (`--min-len` to `--max-small`, default 1 kb-1 Mb) requires **tight local
embedding**: both immediately flanking backbone chains present within
`--max-gap` (100 kb) of the candidate, the candidate's query midpoint lying
between theirs, and its query span within 0.5-2x its reference span. Measured
against the two-human-assembly control, this separates signal from noise far
better than a position-only test -- chimpanzee 96 calls vs control 3 (ratio 32),
where a scaled syntenic-offset test gives 129 vs 44 (ratio 2.9).

*pericentric* (>= `--max-small`) cannot use that test at all: these span tens of
megabases, so no backbone chain lies within 100 kb of both ends, and the parent
chain usually spans the whole event. They are instead required to be
position-consistent -- an inversion flips orientation but preserves position, so
the backbone is interpolated to predict a query midpoint and the candidate must
land within `--offset-tol-frac` x its span (floor `--offset-tol-min`).

Chain nesting is deliberately not resolved by a coverage-based "best chain wins"
assignment: the parent's span covers the gap, so it would hide the inversions.

Outputs one row per merged inverted interval in reference coordinates, plus a
BED of the syntenic backbone used (which regions were assayable at all).
"""
import argparse

import numpy as np
import pandas as pd

import ape_lib as A

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--headers", required=True)
ap.add_argument("--alias", required=True)
ap.add_argument("--ref", default="hg38")
ap.add_argument("--query", required=True)
ap.add_argument("--min-len", type=int, default=1_000,
                help="minimum chain target span (default 1 kb)")
ap.add_argument("--max-small", type=int, default=1_000_000,
                help="boundary between the small and pericentric categories")
ap.add_argument("--max-gap", type=int, default=100_000,
                help="small mode: furthest a flanking backbone chain may sit")
ap.add_argument("--span-ratio", type=float, default=2.0,
                help="small mode: query span must be within [1/r, r] of reference span")
ap.add_argument("--flank-min", type=int, default=10_000,
                help="minimum chain span to count as syntenic backbone")
ap.add_argument("--offset-tol-frac", type=float, default=0.5)
ap.add_argument("--offset-tol-min", type=int, default=100_000)
ap.add_argument("--out", required=True)
ap.add_argument("--backbone-out")
a = ap.parse_args()

alias = pd.read_csv(a.alias, sep="\t")
name = {(r.assembly, r.accession): r.chrom for r in alias.itertuples()}

d = pd.read_csv(a.headers, sep=r"\s+", header=None, names=A.CHAIN_COLS)
d = d[d.tName.isin(A.PRIMARY_CHROMS)].copy()
d["qChrom"] = [name.get((a.query, q), q) for q in d.qName]
d["tLen"] = d.tEnd - d.tStart
# chain files store minus-strand query coordinates in the reverse-complemented frame
d["qF_start"] = np.where(d.qStrand == "-", d.qSize - d.qEnd, d.qStart)
d["qF_end"] = np.where(d.qStrand == "-", d.qSize - d.qStart, d.qEnd)
d["tMid"] = (d.tStart + d.tEnd) // 2
d["qMid"] = (d.qF_start + d.qF_end) // 2

# orthologous query chromosome = the one carrying the most chain score
tot = d.groupby(["tName", "qChrom"]).score.sum().reset_index()
ortho = tot.loc[tot.groupby("tName").score.idxmax()].set_index("tName").qChrom.to_dict()

species = A.SPECIES.get(a.query, a.query)
rows, backbone, n_raw, n_dropped = [], [], 0, 0

for t in A.PRIMARY_CHROMS:
    q = ortho.get(t)
    if q is None:
        continue
    sub = d[(d.tName == t) & (d.qChrom == q)]
    by = sub.groupby("qStrand").score.sum()
    dom = "+" if by.get("+", 0) >= by.get("-", 0) else "-"

    bb = sub[(sub.qStrand == dom) & (sub.tLen >= a.flank_min)].sort_values("tMid")
    for r in bb.itertuples():
        backbone.append((t, r.tStart, r.tEnd, species))

    others = sub[sub.qStrand != dom]
    bb_ts = np.sort(bb.tStart.values)
    bb_qm_by_start = bb.qMid.values[np.argsort(bb.tStart.values)]
    order_end = np.argsort(bb.tEnd.values)
    bb_te = bb.tEnd.values[order_end]
    bb_qm_by_end = bb.qMid.values[order_end]

    picked = []

    # --- small inversions: tight local embedding ---------------------------- #
    small = others[(others.tLen >= a.min_len) & (others.tLen < a.max_small)]
    n_raw += len(small)
    if len(bb) >= 2 and not small.empty:
        qspan = (small.qF_end - small.qF_start).clip(lower=1)
        r = qspan / small.tLen.clip(lower=1)
        small = small[(r >= 1 / a.span_ratio) & (r <= a.span_ratio)]
        for c in small.itertuples():
            li = np.searchsorted(bb_te, c.tStart, "right") - 1
            ri = np.searchsorted(bb_ts, c.tEnd, "left")
            if li < 0 or ri >= bb_ts.size:
                continue
            if c.tStart - bb_te[li] > a.max_gap or bb_ts[ri] - c.tEnd > a.max_gap:
                continue
            lo_q = min(bb_qm_by_end[li], bb_qm_by_start[ri])
            hi_q = max(bb_qm_by_end[li], bb_qm_by_start[ri])
            if lo_q <= c.qMid <= hi_q:
                picked.append((c, "small", float(min(abs(c.qMid - lo_q), abs(c.qMid - hi_q)))))
    n_dropped += len(small) - sum(1 for p in picked if p[1] == "small")

    # --- pericentric inversions: position consistency ------------------------ #
    large = others[others.tLen >= a.max_small]
    n_raw += len(large)
    if len(bb) >= 2 and not large.empty:
        expected = np.interp(large.tMid.values, bb.tMid.values, bb.qMid.values)
        offset = np.abs(large.qMid.values - expected)
        tol = np.maximum(a.offset_tol_frac * large.tLen.values, a.offset_tol_min)
        for c, off, ok in zip(large.itertuples(), offset, offset <= tol):
            if ok:
                picked.append((c, "pericentric", float(off)))
        n_dropped += int((offset > tol).sum())

    if not picked:
        continue

    for category in ("small", "pericentric"):
        sel = sorted((p for p in picked if p[1] == category), key=lambda p: p[0].tStart)
        if not sel:
            continue
        cur = None
        for c, _, off in sel:
            if cur is None:
                cur = [c.tStart, c.tEnd, c.qF_start, c.qF_end, c.score, 1, off]
            elif c.tStart <= cur[1]:
                cur[1] = max(cur[1], c.tEnd)
                cur[2] = min(cur[2], c.qF_start)
                cur[3] = max(cur[3], c.qF_end)
                cur[4] += c.score
                cur[5] += 1
                cur[6] = min(cur[6], off)
            else:
                rows.append((species, a.query, t, *cur[:2], q, *cur[2:], dom, category))
                cur = [c.tStart, c.tEnd, c.qF_start, c.qF_end, c.score, 1, off]
        rows.append((species, a.query, t, *cur[:2], q, *cur[2:], dom, category))


inv = pd.DataFrame(rows, columns=[
    "species", "query_assembly", "ref_chrom", "ref_start", "ref_end", "query_chrom",
    "query_start", "query_end", "chain_score", "n_chains", "colinearity_margin",
    "dominant_orientation", "category"])
inv["ref_span"] = inv.ref_end - inv.ref_start
inv = inv.sort_values(["ref_chrom", "ref_start"], ignore_index=True)

print(f"=== {species} ({a.query}) vs {a.ref} ===")
print(f"candidate chains >= {a.min_len:,} bp: {n_raw:,}   "
      f"dropped by colinearity: {n_dropped:,}")
print(f"merged inversions: {len(inv):,}  ({inv.ref_span.sum()/1e6:.1f} Mb of {a.ref})")
for cat, g in inv.groupby("category"):
    print(f"  {cat:<12} {len(g):>5,}  median {g.ref_span.median()/1e3:,.1f} kb  "
          f"total {g.ref_span.sum()/1e6:.1f} Mb")

inv.to_csv(a.out, sep="\t", index=False)
print(f"wrote {a.out}")

if a.backbone_out:
    bb = pd.DataFrame(backbone, columns=["chrom", "start", "end", "species"])
    bb = bb.sort_values(["chrom", "start"], ignore_index=True)
    bb.to_csv(a.backbone_out, sep="\t", index=False, header=False)
    print(f"wrote {a.backbone_out} ({len(bb):,} backbone chains)")
