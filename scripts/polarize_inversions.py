#!/usr/bin/env python
"""Polarise inversion loci onto branches of the ape tree by Fitch parsimony.

Each per-species table says only whether a segment is inverted *relative to
GRCh38*; that is a pairwise statement, not a direction. Placing the event on a
branch needs the whole pattern across species plus the tree, with siamang as the
outgroup.

Loci are built by merging overlapping calls across species; a species counts as
inverted at a locus when its call and the locus reciprocally overlap by at least
--min-overlap. A species whose syntenic backbone does not cover the locus is
scored missing ('?') rather than not-inverted -- absence of alignment is not
evidence of shared orientation. GRCh38 is state 0 by definition, being the
reference every call was made against.

Fitch's algorithm then reconstructs ancestral states and reports the branches
carrying a change. The number of changes doubles as a quality signal: a real
inversion is usually a single event on one branch, whereas alignment artefacts
give incoherent patterns needing several independent changes. Loci needing more
than one change are the recurrent/homoplastic class -- biologically the
interesting "toggling" regions, but also where artefacts accumulate.
"""
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd

import ape_lib as A

# Guide tree of the alignment (README), siamang as outgroup.
TREE = ('root', [
    ('Hominidae', [
        ('Pongo', ['Sumatran_orangutan', 'Bornean_orangutan']),
        ('African_apes', [
            'gorilla',
            ('Homininae', [
                ('Pan', ['bonobo', 'chimpanzee']),
                ('Homo', ['human_CHM13', 'human_GRCh38']),
            ]),
        ]),
    ]),
    'siamang',
])
TAXA = ['Sumatran_orangutan', 'Bornean_orangutan', 'gorilla', 'bonobo',
        'chimpanzee', 'human_CHM13', 'human_GRCh38', 'siamang']


def children(node):
    return node[1] if isinstance(node, tuple) else []


def label(node):
    return node[0] if isinstance(node, tuple) else node


def fitch(node, states):
    """Down-pass: returns (state_set, n_changes_below)."""
    if not isinstance(node, tuple):
        return states[node], 0
    sets, changes = [], 0
    for c in children(node):
        s, ch = fitch(c, states)
        sets.append(s)
        changes += ch
    inter = set.intersection(*sets)
    if inter:
        return inter, changes
    return set.union(*sets), changes + 1


def assign(node, states, parent_state, out):
    """Up-pass: fix states and record the branches carrying a change."""
    s, _ = fitch(node, states)
    me = parent_state if parent_state in s else min(s)
    if parent_state is not None and me != parent_state:
        out.append(label(node))
    for c in children(node):
        assign(c, states, me, out)
    return me


def merge_intervals(starts, ends):
    """Union of a set of intervals, returned sorted and non-overlapping.

    Backbone chains nest and overlap heavily, so they must be unioned before
    they can be searched: after sorting by start, `ends` is NOT monotonic, and
    a searchsorted over it silently returns nonsense.
    """
    if starts.size == 0:
        return starts, ends
    o = np.argsort(starts, kind='stable')
    s, e = starts[o], ends[o]
    out_s, out_e = [], []
    cs, ce = s[0], e[0]
    for i in range(1, s.size):
        if s[i] <= ce:
            ce = max(ce, e[i])
        else:
            out_s.append(cs); out_e.append(ce)
            cs, ce = s[i], e[i]
    out_s.append(cs); out_e.append(ce)
    return np.asarray(out_s), np.asarray(out_e)


def coverage_fraction(iv_start, iv_end, starts, ends):
    """Fraction of [iv_start, iv_end) covered by merged, sorted intervals."""
    if starts.size == 0:
        return 0.0
    i = np.searchsorted(ends, iv_start, 'right')
    j = np.searchsorted(starts, iv_end, 'left')
    if i >= j:
        return 0.0
    cov = np.minimum(ends[i:j], iv_end) - np.maximum(starts[i:j], iv_start)
    return float(cov[cov > 0].sum()) / max(iv_end - iv_start, 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", nargs="+", required=True)
    ap.add_argument("--backbones", nargs="+", required=True)
    ap.add_argument("--min-overlap", type=float, default=0.5)
    ap.add_argument("--min-backbone", type=float, default=0.5,
                    help="backbone coverage below this scores the species as missing")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    calls = pd.concat([pd.read_csv(p, sep="\t") for p in a.tables], ignore_index=True)
    calls = calls[calls.ref_span > 0]
    if "category" not in calls.columns:
        calls["category"] = "small"
    print(f"input calls: {len(calls):,} across {calls.species.nunique()} species")
    print(calls.category.value_counts().to_string())

    bb = {}
    for p in a.backbones:
        b = pd.read_csv(p, sep="\t", header=None,
                        names=["chrom", "start", "end", "species"])
        sp = b.species.iloc[0]
        bb[sp] = {c: merge_intervals(g.start.values, g.end.values)
                  for c, g in b.groupby("chrom")}

    # loci = union-merge of all calls, built separately per category so that a
    # pericentric event never swallows the small inversions inside it
    loci = []
    for category, cat_calls in calls.groupby("category"):
      for chrom, g in cat_calls.sort_values(["ref_chrom", "ref_start"]).groupby("ref_chrom"):
        cs = ce = None
        for r in g.itertuples():
            if cs is None:
                cs, ce = r.ref_start, r.ref_end
            elif r.ref_start <= ce:
                ce = max(ce, r.ref_end)
            else:
                loci.append((chrom, cs, ce, category)); cs, ce = r.ref_start, r.ref_end
        loci.append((chrom, cs, ce, category))
    print(f"merged loci: {len(loci):,}")

    by_species = {(sp, cat): {c: g[["ref_start", "ref_end"]].values
                              for c, g in gg.groupby("ref_chrom")}
                  for (sp, cat), gg in calls.groupby(["species", "category"])}

    rows = []
    for chrom, lo, hi, category in loci:
        L = hi - lo
        states, cov = {}, {}
        for sp in TAXA:
            if sp == 'human_GRCh38':
                states[sp], cov[sp] = {0}, 1.0
                continue
            iv = by_species.get((sp, category), {}).get(chrom)
            inverted = False
            if iv is not None:
                ov = np.minimum(iv[:, 1], hi) - np.maximum(iv[:, 0], lo)
                ov = np.maximum(ov, 0)
                call_len = iv[:, 1] - iv[:, 0]
                inverted = bool((((ov / L) >= a.min_overlap) &
                                 ((ov / np.maximum(call_len, 1)) >= a.min_overlap)).any())
            s, e = bb.get(sp, {}).get(chrom, (np.array([]), np.array([])))
            f = coverage_fraction(lo, hi, s, e)
            cov[sp] = f
            if inverted:
                states[sp] = {1}
            elif f >= a.min_backbone:
                states[sp] = {0}
            else:
                states[sp] = {0, 1}          # not assayable -> missing

        branches = []
        assign(TREE, states, None, branches)
        _, n_changes = fitch(TREE, states)
        n_missing = sum(1 for sp in TAXA if states[sp] == {0, 1})
        rows.append((chrom, lo, hi, L, category,
                     *[('?' if states[sp] == {0, 1} else str(min(states[sp]))) for sp in TAXA],
                     n_changes, ";".join(branches) if branches else "none", n_missing,
                     round(min(cov.values()), 3)))

    out = pd.DataFrame(rows, columns=["chrom", "start", "end", "span", "category", *TAXA,
                                      "n_changes", "branches", "n_missing",
                                      "min_backbone_cov"])
    out = out.sort_values(["chrom", "start"], ignore_index=True)

    for cat, g in out.groupby("category"):
        clean = g[g.n_changes == 1]
        print(f"\n=== {cat}: {len(g):,} loci, {len(clean):,} resolve to one branch ===")
        print(clean.branches.value_counts().to_string())

    out.to_csv(a.out, sep="\t", index=False)
    print(f"\nwrote {a.out}")


if __name__ == '__main__':
    main()
