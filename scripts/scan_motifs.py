#!/usr/bin/env python
"""Call the CTCF motif orientation of every site on one hg38 chromosome (README.md §9.2).

Scores the JASPAR MA0139.1 PWM at every position of the chromosome on both
strands -- the same scoring as the chr21 validation in §9 -- and reads off, per
site, the best hit within +/-100 bp of the summit and the number of hits on each
strand across the whole site interval. The output is keyed by `chrom`/`start`
and merged back into the site table by `annotate_motifs.py`.
"""
import argparse
import json

import pandas as pd
import ctcf_lib as L

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--sites", required=True, help="consensus site table (.tsv.gz)")
ap.add_argument("--fasta", required=True, help="hg38 chromosome .fa.gz")
ap.add_argument("--pfm", required=True, help="JASPAR MA0139.1 JSON")
ap.add_argument("--chrom", required=True)
ap.add_argument("--flank", type=int, default=L.MOTIF_FLANK)
ap.add_argument("--cutoff", type=float, default=L.MOTIF_REL_CUTOFF)
ap.add_argument("--out", required=True)
a = ap.parse_args()

S = L.encode_sequence(L.read_fasta_gz(a.fasta))
lo, rc, mn, mx = L.load_pwm(json.load(open(a.pfm))["pfm"])
W = lo.shape[1]
print(f"{a.chrom}: {S.size:,} bp, motif width {W}")

fwd, rev = L.scan_pwm(S, lo, rc, mn, mx)
print(f"  windows >= {a.cutoff}: {(fwd >= a.cutoff).sum():,} plus, {(rev >= a.cutoff).sum():,} minus")

sites = pd.read_csv(a.sites, sep="\t", usecols=["chrom", "start", "end", "summit"])
sub = sites[sites.chrom == a.chrom]
ann = L.annotate_motifs(sub, fwd, rev, W, a.flank, a.cutoff)
out = pd.concat([sub[["chrom", "start"]], ann], axis=1)

n = len(out)
plus, minus = (out.motif_strand == "+").sum(), (out.motif_strand == "-").sum()
print(f"  {n:,} sites: {plus:,} oriented +, {minus:,} oriented -, {n - plus - minus:,} no motif")
print(f"  motifs in both orientations: {((out.n_motifs_plus > 0) & (out.n_motifs_minus > 0)).sum():,}")

# explicit, so the write-then-move scratch name cannot change the format
out.to_csv(a.out, sep="\t", index=False, compression="gzip")
print(f"wrote {a.out}")
