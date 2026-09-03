# Scripts

Steps invoked by the GWF workflow in `../workflow.py`. Each is a standalone CLI
(`--help` works); `ctcf_lib.py` holds the shared logic and is imported, not run.
The specification they implement is `../METHODS.md` — the section numbers in the
docstrings point there.

| script | stage |
|---|---|
| `ctcf_lib.py` | shared: ENCODE selection rules, narrowPeak parsing, clustering, chain liftover, motif scoring |
| `fetch_encode_metadata.py` | retrieval — ENCODE batch-metadata table (§2.1) |
| `select_encode_files.py` | retrieval — one peak file per experiment → provenance manifest (§2.2–2.3) |
| `fetch_peak_chunk.py` | retrieval — download and parse one chunk of peak files (§3) |
| `fetch_reference.py` | retrieval — liftover chains, hg38 chr21, JASPAR MA0139.1, CHM13 chrom sizes |
| `build_consensus.py` | compilation — cluster pooled summits into consensus sites (§4–5) |
| `liftover_chm13.py` | compilation — append T2T-CHM13v2.0 coordinates, keeping strand (§11) |
| `export_parquet.py` | publication — shard a site table into a parquet dataset of sub-40 MB part files (§6.1) |
| `ape_lib.py` | shared: T2T ape species map, guide-tree order, chain/chromAlias URLs |
| `fetch_ape_chain_headers.py` | retrieval — stream one all-to-all chain file, keep header lines (§12.1) |
| `fetch_chromalias.py` | retrieval — GenBank accession → chromosome name maps |
| `call_inversions.py` | compilation — call inversions against the dominant orientation per pair (§12.2) |
| `polarize_inversions.py` | compilation — place inversion loci on tree branches by Fitch parsimony (§13) |
| `merge_inversions.py` | compilation — combine the per-species inversion tables (§12.3) |
| `validate_motif.py` | validation — CTCF motif enrichment on chr21 vs shuffled controls (§9) |
| `validate_liftover.py` | validation — CHM13 coordinate checks and orientation flips (§11.5–11.6) |
| `dump_summits.py` + `merge_sites.awk` + `compare_consensus.py` | validation — re-derive one chromosome with an independent awk implementation and diff it (§7) |

The validation scripts exit non-zero when a check fails. The workflow writes
their report to a scratch path and only moves it into `results/` on success, so
a failed check leaves the target incomplete instead of publishing a bad report —
read the details with `gwf logs validate_motif`.

The other files here (`rename.py`, `git-commit.sh`, `docs-run-notebooks.sh`,
`test.sh`) are project scaffolding, unrelated to the pipeline.
