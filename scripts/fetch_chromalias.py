#!/usr/bin/env python
"""Fetch UCSC chromAlias tables mapping GenBank accessions to chromosome names.

The chain files name ape chromosomes by GenBank accession (CM054434.2); the
alias table turns those into the assembly's own chromosome names.
"""
import argparse

import requests

import ape_lib as A

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--out", required=True)
a = ap.parse_args()

rows = []
for acc in A.GCA_ASSEMBLIES:
    r = requests.get(A.chromalias_url(acc), timeout=120)
    r.raise_for_status()
    for line in r.text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        rows.append(f"{acc}\t{f[0]}\t{f[-1]}")
    print(f"{acc}: {len(rows)} rows so far")

with open(a.out, "w") as fh:
    fh.write("assembly\taccession\tchrom\n")
    fh.write("\n".join(rows) + "\n")
print(f"wrote {len(rows)} aliases to {a.out}")
