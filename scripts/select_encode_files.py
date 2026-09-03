#!/usr/bin/env python
"""Pick one peak file per ENCODE experiment (README.md §2.2-2.3).

Newest released analysis, best IDR output type, largest file as tie-break;
experiments whose selected file carries a read-depth ERROR audit are dropped
unless --include-low-quality is given.
"""
import argparse
import pandas as pd
import ctcf_lib as L

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--metadata", required=True)
ap.add_argument("--out", required=True, help="destination manifest TSV")
ap.add_argument("--include-low-quality", action="store_true")
a = ap.parse_args()

meta = pd.read_csv(a.metadata, sep="\t", low_memory=False)
sel = L.select_files(meta, a.include_low_quality)

print(f"metadata rows            {len(meta):,} "
      f"({meta['Experiment accession'].nunique()} experiments)")
print(f"dropped for ERROR audits {sel.attrs['n_dropped']}")
print(f"selected                 {len(sel)} experiments, "
      f"{sel['biosample'].nunique()} biosamples")
print(f"download volume          {sel['size_bytes'].sum()/1e6:.0f} MB")

sel.to_csv(a.out, sep="\t", index=False)
print(f"wrote {a.out}")
