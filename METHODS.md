# ENCODE CTCF binding sites, hg38 — compilation methods

**Dataset:** `ctcf_sites_hg38.tsv.gz` — 337,104 consensus CTCF binding sites on the
GRCh38/hg38 primary assembly, derived from every released human CTCF ChIP-seq
experiment in ENCODE.

**ENCODE snapshot:** 2026-09-02 (the portal is a living database; see §10.7).

This document specifies the compilation completely enough to reproduce the dataset
from scratch without the accompanying Python code. §7 gives a verified shell-only
recipe; §8 gives the numbers a reimplementation must reproduce.

---

## 0. Where things live

This document specifies the pipeline implemented by `workflow.py` and the
scripts in `scripts/`. Section numbers referenced from the code point here.

| path | what it is |
|---|---|
| `workflow.py` | the GWF workflow: retrieval, compilation, validation |
| `scripts/` | the steps the workflow invokes (see `scripts/README.md`) |
| `results/ctcf_sites_hg38/` | **the dataset**: 337,104 CTCF sites, parquet (§6.1) |
| `results/ctcf_sites_hg38_chm13/` | the same sites with CHM13 coordinates added, parquet (§11) |
| `results/ctcf_encode_files.tsv` | provenance manifest: the 435 ENCODE files used, with md5sums (§6.2) |
| `results/validation/` | the three validation reports (§7, §9, §11.5) |
| `steps/` | intermediates — metadata, peak chunks, reference files, the TSV site tables (git-ignored) |

Run the whole thing with:

```bash
pixi run gwf workers -n 4 &      # local backend only
pixi run gwf run
```

Loading the dataset needs nothing but pandas:

```python
from pd_lfs import read_parquet
sites = read_parquet("results/ctcf_sites_hg38")
sites = sites[sites.n_experiments >= 10]      # see §10.1 before choosing a threshold
```

`scripts/ctcf_lib.py` exposes `load_ctcf_sites(path, min_experiments=...)`, which
accepts a parquet dataset (a directory or an https URL to one) or a `.tsv.gz`,
and handles both the hg38-only and the dual-coordinate table.

---

## 1. Summary of the compilation

| Stage | Result |
|---|---|
| ENCODE files matching the query | 1,971 files across 457 experiments |
| Candidate peak files after quality filters | 1,066 files across 457 experiments |
| One best peak file per experiment | 457 files |
| After excluding read-depth ERROR audits | **435 files / 435 experiments / 158 biosamples** |
| Peak summits pooled | 18,521,082 |
| Consensus sites after clustering | **337,104** |
| Genome covered | 93.0 Mb (3.00% of hg38) |

---

## 2. Source data

### 2.1 Enumerating candidate files

All released human CTCF TF ChIP-seq experiments with GRCh38 narrowPeak output were
enumerated through the ENCODE portal's batch-metadata endpoint:

```
https://www.encodeproject.org/metadata/?type=Experiment&assay_title=TF+ChIP-seq&target.label=CTCF&assembly=GRCh38&status=released&files.file_type=bed+narrowPeak&files.assembly=GRCh38
```

This returns a tab-separated table with one row per file (1,971 rows, 457 distinct
experiment accessions on the snapshot date). The columns used downstream are:

`File accession`, `Output type`, `File assembly`, `Experiment accession`,
`Biosample term name`, `Biosample type`, `Biosample organism`,
`Experiment date released`, `Size`, `md5sum`, `File download URL`, `File Status`,
`File analysis title`, `File analysis status`, `Audit ERROR`, `Audit NOT_COMPLIANT`.

### 2.2 Selecting one peak file per experiment

An experiment can carry several peak files (different pipeline versions, different IDR
flavours). Rows were first restricted to:

- `Output type` ∈ {`IDR thresholded peaks`, `optimal IDR thresholded peaks`,
  `conservative IDR thresholded peaks`, `pseudoreplicated IDR thresholded peaks`}
  — this excludes `peaks and background as input for IDR` and unthresholded `peaks`,
  which are not statistically thresholded peak calls;
- `File Status` == `released`;
- `File analysis status` == `released` (excludes superseded/archived analyses);
- `File assembly` == `GRCh38`;
- `Biosample organism` == `Homo sapiens`.

1,066 files survive, still covering all 457 experiments. Within each experiment the
remaining files were ranked and the top one taken, sorting by:

1. **Analysis version**, ascending rank — newest pipeline first:

   | rank | `File analysis title` |
   |---|---|
   | 0 | ENCODE4 v1.8.0 GRCh38 |
   | 1 | ENCODE4 v1.7.0 GRCh38 |
   | 2 | ENCODE4 v1.6.1 GRCh38 |
   | 3 | ENCODE4 v1.5.1 GRCh38 |
   | 4 | ENCODE4 v1.5.0 GRCh38 |
   | 5 | ENCODE4 v1.4.0 GRCh38 |
   | 6 | ENCODE3 GRCh38 |
   | 7 | Lab custom GRCh38 |
   | 99 | anything else / missing |

2. **Output type**, ascending rank: `IDR thresholded peaks` (0) <
   `optimal IDR thresholded peaks` (1) < `conservative IDR thresholded peaks` (2) <
   `pseudoreplicated IDR thresholded peaks` (3).

3. **`Size`**, descending, as the final tie-break.

### 2.3 Quality exclusion

An experiment was dropped if the selected file's `Audit ERROR` field contained
(case-insensitive substring match) either `extremely low read depth` or
`control extremely low read depth`. These indicate the peak set itself is
untrustworthy. **22 experiments were dropped**, leaving **435**.

Other ENCODE audits were *retained but recorded*, not used to exclude: `missing
control alignments` is a metadata gap rather than a data-quality problem, and
NOT_COMPLIANT flags (`insufficient read depth`, `severe bottlenecking`, `poor library
complexity`, `partially characterized antibody`, `unreplicated experiment`, …) mark
sub-optimal but usable experiments. Both audit strings are carried verbatim in the
`audit_error` / `audit_not_compliant` columns of `ctcf_encode_files.tsv`, so any
stricter exclusion can be applied post hoc.

### 2.4 Composition of the selected set

All 435 selected files turned out to be `IDR thresholded peaks` from ENCODE4
analyses — every experiment, including the older ENCODE3 submissions, has been
reprocessed with the ENCODE4 uniform pipeline. So the catalogue is uniformly
processed; no ENCODE3 or lab-custom peak calls contribute.

| Analysis | files | | Biosample type | files |
|---|---|---|---|---|
| ENCODE4 v1.5.1 GRCh38 | 267 | | tissue | 244 |
| ENCODE4 v1.8.0 GRCh38 | 126 | | cell line | 122 |
| ENCODE4 v1.5.0 GRCh38 | 31 | | primary cell | 42 |
| ENCODE4 v1.6.1 GRCh38 | 7 | | in vitro differentiated cells | 20 |
| ENCODE4 v1.7.0 GRCh38 | 3 | | organoid | 7 |
| ENCODE4 v1.4.0 GRCh38 | 1 | | | |

158 distinct biosamples. The most heavily replicated are dorsolateral prefrontal
cortex (59 experiments), HCT116 (19), A549 (18), heart left ventricle (14), spleen
(13); see §10.3.

The exact file list — experiment accession, file accession, biosample, pipeline
version, size, **md5sum** and download URL — is shipped as `ctcf_encode_files.tsv`
(435 rows). Total download: 310 MB.

---

## 3. From peak calls to summits

Each selected file is a gzipped ENCODE **narrowPeak** (BED6+4), 0-based half-open
coordinates, tab-separated, no header:

| col | field | used |
|---|---|---|
| 1 | `chrom` | yes |
| 2 | `chromStart` (0-based) | yes |
| 3 | `chromEnd` (exclusive) | yes |
| 4 | `name` | no |
| 5 | `score` | no |
| 6 | `strand` | no |
| 7 | `signalValue` — fold enrichment over control | yes |
| 8 | `-log10(pValue)` | no |
| 9 | `-log10(qValue)` | yes |
| 10 | `peak` — summit offset from `chromStart`, `-1` if not called | yes |

For every peak a single **summit** position was computed:

```
summit = chromStart + peak            if peak >= 0
summit = floor((chromStart + chromEnd) / 2)   otherwise
```

Peaks were restricted to the **primary assembly**: `chr1`–`chr22`, `chrX`, `chrY`.
This excludes `chrM`, unplaced/unlocalized scaffolds and alt contigs. In practice it
is a no-op safeguard — ENCODE4 IDR peak files are already confined to these 24
sequences (verified: 0 of 510,660 peaks in a 12-file sample fell outside them).

Each retained summit carries four attributes into the next stage: the **experiment
accession** it came from, the **biosample term name** of that experiment, its
`signalValue`, and its `-log10(qValue)`.

Pooling all 435 files gives **18,521,082 summits** (mean 42,577 peaks per
experiment).

---

## 4. Clustering summits into consensus sites

Independent experiments call the same physical CTCF site at slightly offset
positions, so summits were clustered across experiments.

**Sort.** All pooled summits sorted by chromosome, then by summit position ascending.
(Chromosome order only has to be *consistent*; it does not affect which summits group
together.)

**Cluster rule — single linkage at 150 bp.** Walking the sorted list, a new cluster
starts at row *i* when

```
chrom[i] != chrom[i-1]   OR   summit[i] - summit[i-1] > 150
```

i.e. a maximal run of summits in which every consecutive pair is ≤ 150 bp apart forms
one site. Equivalently: place a 151 bp window `[summit-75, summit+76)` around every
summit and take connected components of overlapping windows. `150 = 2 × slop` with
`slop = 75`; the parameter is `--slop` in the accompanying script.

**Site interval.** For a cluster spanning summits `s_min … s_max`:

```
start = max(s_min - 75, 0)
end   = s_max + 75 + 1
```

0-based, half-open, so `width = end - start`, and a cluster containing a single summit
yields a 151 bp site. No site exceeds its chromosome length in hg38 (verified against
`hg38.chrom.sizes`), so no right-clipping was applied; `start` is clipped at 0.

**Rationale and chaining.** 150 bp is generous relative to the cross-experiment
summit jitter of a sharp factor like CTCF but tight enough that single-linkage
chaining stays local. Observed widths: minimum 151 bp, median 182 bp, mean 276 bp,
90th percentile 514 bp, 99th percentile 938 bp, maximum 2,719 bp; only 2,641 sites
(0.78%) exceed 1 kb and 43 exceed 2 kb. Those long sites are genuinely clustered CTCF
arrays, but if strictly one-motif-per-site elements are required, either filter on
`width` or rebuild with a smaller `--slop`.

---

## 5. Per-site annotation

For each cluster, over the pooled peaks belonging to it:

| column | definition |
|---|---|
| `chrom` | chromosome |
| `start` | `max(min(summit) - 75, 0)`, 0-based inclusive |
| `end` | `max(summit) + 76`, exclusive |
| `summit` | summit of the peak with the **highest** `signalValue` in the cluster; ties broken by lowest genomic coordinate |
| `n_experiments` | number of **distinct experiment accessions** contributing ≥ 1 peak |
| `n_biosamples` | number of **distinct biosample term names** contributing ≥ 1 peak |
| `frac_experiments` | `n_experiments / 435` |
| `max_signal` | maximum `signalValue` (narrowPeak col 7) |
| `mean_signal` | mean `signalValue` over **all peaks** in the cluster (not per experiment; see §10.2) |
| `max_neglog10q` | maximum `-log10(qValue)` (narrowPeak col 9) |
| `n_peaks` | number of pooled peaks in the cluster (≥ `n_experiments`) |
| `width` | `end - start` |

Rows are sorted by chromosome (`chr1`…`chr22`, `chrX`, `chrY`) then `start`.

---

## 6. Output file specification

### 6.1 `results/ctcf_sites_hg38/`

A parquet dataset written by `pd_lfs.write_parquet`: `part-*.parquet` shards each
below 40 MB — comfortably under GitHub's 50 MB warning — plus a `_manifest.json`
index that lets `pd_lfs.read_parquet` read the dataset over plain HTTPS without
directory listing. 337,104 rows, 12 columns in the order listed in §5. Conversion
is lossless; the column dtypes below are carried in the manifest and restored on
read. The pipeline also emits the same table as a gzipped TSV at
`steps/consensus/ctcf_sites_hg38.tsv.gz`, which is what the shell recipe in §7
reproduces. Coordinates are **0-based half-open (BED convention) on GRCh38/hg38**. Column
types: `chrom` category; `start`/`end`/`summit`/`n_peaks`/`width` int32;
`n_experiments`/`n_biosamples` int16; `frac_experiments`/`max_signal`/`mean_signal`/
`max_neglog10q` float32.

To load as a pandas DataFrame with no other dependency:

```python
import pandas as pd
sites = pd.read_csv("ctcf_sites_hg38.tsv.gz", sep="\t")
```

Converting to 1-based inclusive (GTF/VCF-style) coordinates: `start + 1` … `end`.

### 6.2 `results/ctcf_encode_files.tsv`

Provenance manifest, one header line, 435 rows: `experiment`, `file_accession`,
`biosample`, `biosample_type`, `output_type`, `analysis`, `lab`, `date_released`,
`size_bytes`, `md5sum`, `url`, `audit_error`, `audit_not_compliant`. Every input file
is identified by ENCODE accession and md5, so an independent rebuild can confirm it
used byte-identical inputs.

---

## 7. Reproducing the dataset without the Python code

Requires only `curl`, `gzip`, `sort` and `awk`. This recipe was executed on chr21 and
produces output **identical** to the shipped file: 4,407/4,407 sites with zero
mismatches in `start`, `end`, `summit`, `n_experiments`, `n_biosamples` and `n_peaks`,
and agreement to float32 precision (≤ 1.3e-4) on the signal columns.

Set `LC_ALL=C` throughout: a locale with comma decimal separators corrupts `awk`'s
numeric output, and locale collation changes `sort` order.

**Step 1 — get the file list.** Either use the shipped `ctcf_encode_files.tsv`, or
rebuild it by downloading the metadata TSV of §2.1 and applying §2.2–2.3.

**Step 2 — download the peak files** (310 MB; verify against column 10, `md5sum`):

```bash
mkdir -p peaks
tail -n +2 ctcf_encode_files.tsv | cut -f2,11 | while IFS=$'\t' read -r acc url; do
    [ -s "peaks/$acc.bed.gz" ] || curl -sL -o "peaks/$acc.bed.gz" "$url"
done
```

**Step 3 — emit annotated summits** (§3):

```bash
tail -n +2 ctcf_encode_files.tsv | cut -f1,2,3 | while IFS=$'\t' read -r expt acc bios; do
    gzip -dc "peaks/$acc.bed.gz" | LC_ALL=C awk -v e="$expt" -v b="$bios" '
        BEGIN { FS = OFS = "\t" }
        $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {
            s = ($10 >= 0) ? $2 + $10 : int(($2 + $3) / 2)
            printf "%s\t%d\t%d\t%s\t%s\t%.6f\t%.6f\n", $1, s, s + 1, e, b, $7, $9
        }'
done > summits.bed                       # 18,521,082 lines
```

**Step 4 — sort and cluster** (§4–5), using the shipped `scripts/merge_sites.awk`:

```bash
LC_ALL=C sort -k1,1 -k2,2n summits.bed > summits.sorted.bed
LC_ALL=C awk -f merge_sites.awk summits.sorted.bed > sites.raw.tsv
```

`scripts/merge_sites.awk` is a single streaming pass implementing §4 and §5 exactly; its
`SLOP` variable (default 75) can be overridden with `awk -v SLOP=50`. It emits
`chrom, start, end, summit, n_experiments, n_biosamples, n_peaks, max_signal,
mean_signal, max_neglog10q`.

**Step 5 — add the derived columns and the header:**

```bash
{ printf 'chrom\tstart\tend\tsummit\tn_experiments\tn_biosamples\tfrac_experiments\tmax_signal\tmean_signal\tmax_neglog10q\tn_peaks\twidth\n'
  LC_ALL=C awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4,$5,$6,$5/435,$8,$9,$10,$7,$3-$2}' sites.raw.tsv
} | gzip > ctcf_sites_hg38.tsv.gz
```

Note this yields rows in lexicographic chromosome order (`chr1, chr10, chr11, …`)
whereas the shipped file uses karyotype order (`chr1, chr2, …, chr22, chrX, chrY`);
the row *set* is identical.

If you prefer bedtools, `bedtools merge -d 149 -c 4,5,7,6,6,7 -o
count_distinct,count_distinct,count,max,mean,max` on the sorted 1 bp summit BED,
followed by `start -= 75; end += 75`, gives the same intervals and counts — `-d 149`
because bedtools measures the gap between book-ended features as 0. It cannot
reproduce the `summit` column (arg-max of `signalValue`), so `merge_sites.awk` is the
faithful route.

---

## 8. Checkpoints for a reimplementation

Any faithful rebuild from the same ENCODE snapshot must reproduce these exactly.

**Totals:** 435 input experiments; 18,521,082 pooled summits; 337,104 sites;
93.0 Mb covered; `sum(start) = 26,278,219,618,224`; `sum(n_experiments) = 16,250,508`.

**Sites per chromosome:**

| chr | sites | chr | sites | chr | sites | chr | sites |
|---|---|---|---|---|---|---|---|
| chr1 | 30,478 | chr7 | 17,936 | chr13 | 8,588 | chr19 | 11,249 |
| chr2 | 26,462 | chr8 | 15,773 | chr14 | 10,899 | chr20 | 12,776 |
| chr3 | 21,340 | chr9 | 14,517 | chr15 | 11,105 | chr21 | 4,407 |
| chr4 | 15,759 | chr10 | 15,793 | chr16 | 11,769 | chr22 | 6,664 |
| chr5 | 17,999 | chr11 | 16,508 | chr17 | 16,491 | chrX | 7,938 |
| chr6 | 18,827 | chr12 | 16,011 | chr18 | 7,420 | chrY | 395 |

**Recurrence distribution** (`n_experiments`; mean 48.2, median 2, max 431):

| threshold | ≥1 | ≥2 | ≥5 | ≥10 | ≥50 | ≥100 | ≥200 | ≥300 | ≥400 |
|---|---|---|---|---|---|---|---|---|---|
| sites | 337,104 | 194,177 | 135,020 | 107,500 | 61,517 | 46,049 | 32,916 | 24,899 | 13,892 |

142,927 sites (42.4%) are seen in exactly one experiment.

**Signal tracks recurrence** — mean `max_signal` by recurrence bin: 1 experiment 23.6;
2–5 45.9; 6–25 91.4; 26–100 187.6; 101–250 407.3; 251–435 1169.2.

**Widths:** min 151, median 182, mean 275.9, p90 514, p99 938, p99.9 1,485, max 2,719.

---

## 9. Validation

Independent confirmation that the coordinates are hg38 and that the sites are real
CTCF sites, using sequence evidence that never entered the compilation.

**Procedure.** hg38 `chr21` was downloaded from UCSC
(`https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz`) and the
JASPAR CTCF matrix MA0139.1 (19 bp) from
`https://jaspar.elixir.no/api/v1/matrix/MA0139.1/?format=json`. The PFM was converted
to a log2-odds PWM against a uniform 0.25 background with a pseudocount of 1% of the
mean column total, distributed as `0.25 × pseudo` per base. Every 19-mer within
±100 bp of each site's `summit` was scored on both strands (windows containing `N`
were rejected), and each site was scored by its best hit expressed as a JASPAR
*relative score* `(score − min) / (max − min)`, where min/max are the sums of the
per-position minima/maxima of the PWM. A site "has a motif" if that relative score
≥ 0.80. Backgrounds: 5,000 random non-`N` chr21 positions, and the same summits
shifted 5 kb downstream.

**Result:**

| set | motif rate |
|---|---|
| random chr21 positions (n = 5,000) | 7.3% |
| sites in 1 experiment (n = 1,988) | 18.7% |
| sites in 2–9 experiments (n = 1,177) | 40.8% |
| sites in ≥10 experiments (n = 1,242) | 72.5% |
| sites in ≥100 experiments (n = 464) | 88.1% |
| sites in ≥300 experiments (n = 230) | 93.9% |
| ≥100-experiment sites shifted +5 kb (n = 464) | 10.1% |

The monotone rise with recurrence, and the collapse to near-background under a 5 kb
shift, confirm both the assembly and the biological identity of the sites. Wrong-
assembly coordinates would sit at the ~7% background rate throughout. The workflow target
`validate_motif` reruns this end to end on every build.

---

## 10. Interpretation and caveats

**10.1 The table is a union, and `n_experiments` is the column that matters.**
No occupancy threshold was imposed. 42% of sites come from a single experiment, and
§9 shows those are only ~2.5-fold motif-enriched over background — largely noise and
genuinely rare cell-type-specific binding, mixed. Filter to `n_experiments >= 10` for
general use (107,500 sites, 72.5% motif-positive) or `>= 100` for the constitutive
core (46,049 sites, 88.1% motif-positive — consistent with the long-standing estimate
of ~40–50k cell-type-invariant CTCF sites). Keep the singletons only when
specifically studying cell-type-specific binding, and treat them with suspicion.

**10.2 `signalValue` is not comparable across experiments.** ENCODE `signalValue` is
fold enrichment over that experiment's own control, uncorrected for depth, antibody
lot or pipeline version. `max_signal` and `mean_signal` are therefore useful for
*ranking sites within this table* but should not be read as quantitative occupancy,
and `mean_signal` is additionally biased by which experiments happened to call the
site. `n_experiments` is the more robust strength proxy.

**10.3 Biosample representation is uneven.** 435 experiments cover only 158
biosamples; dorsolateral prefrontal cortex alone contributes 59. `n_experiments` thus
partly reflects sampling effort. `n_biosamples` (max 158) is the less biased
recurrence measure — use it when the question is "how many cell types", and
`n_experiments` when the question is "how reproducibly called".

**10.4 No motif orientation.** CTCF loop anchors are directional and the convergent-
motif rule governs loop formation, but this table carries no motif strand. Adding it
requires a PWM scan of the full hg38 sequence (§9 does exactly this for chr21; scaling
it genome-wide needs the ~950 MB hg38 FASTA). Sites in long clusters (§4) may contain
several motifs in different orientations.

**10.5 IDR peak sets are conservative by construction.** IDR thresholding demands
reproducibility between replicates, so weak-but-real sites are missed, particularly in
low-depth experiments. The catalogue is a high-precision, not a high-recall, view of
CTCF binding.

**10.6 chrY is underrepresented** (395 sites): many ENCODE biosamples are female, and
chrY is repeat-rich and poorly mappable. chrX sites (7,938) come from a mixture of
active and inactive X backgrounds.

**10.7 The ENCODE snapshot moves.** The portal continuously adds experiments and
reprocesses old ones with newer pipeline versions, and `File analysis status` flips
from `released` to `archived` as it does. Rebuilding on a later date will select
different files and yield different totals. The dataset shipped here corresponds to
the 2026-09-02 snapshot, pinned exactly by the accessions and md5sums in
`ctcf_encode_files.tsv`.

**10.8 Excluded experiments.** The 22 experiments dropped in §2.3 are not listed in
the manifest. To recover them, repeat §2.2 without the audit filter (`scripts/select_encode_files.py --include-low-quality`) and diff the experiment
accessions.

---

## 11. T2T-CHM13v2.0 coordinates (derived)

`results/ctcf_sites_hg38_chm13/` is the §6.1 table with six columns appended. The
hg38 columns are unchanged and remain authoritative; the CHM13 columns are
**derived by liftover, not by native alignment**.

### 11.1 Why derived

ENCODE holds no CTCF ChIP-seq on CHM13. Queried on 2026-09-02, the assembly
facet for CTCF TF ChIP-seq offers only GRCh38 (457 experiments), hg19 (270),
mm10, dm3 and dm6. Across the *entire* portal (1,659,000 files) the assembly
`T2T-CHM13` covers 48 files, all of them long-read RNA-seq modification calls
from 4 experiments (ENCSR872GND, ENCSR111GJE, ENCSR697CSS, ENCSR543NWW). There
is no ChIP-seq of any target on CHM13.

The only CHM13-native ENCODE CTCF peaks found anywhere are 6 files in the T2T
ENCODE reanalysis (`s3://human-pangenomics/T2T/CHM13/assemblies/annotation/
regulation/ENCODE/macs2_peak/`), covering 22Rv1, C4-2B, RWPE1, RWPE2, VCaP and
prostate epithelium — a prostate-focused reanalysis of 44 peak files total, not
an ENCODE-wide remap. A CHM13-native equivalent of this 435-experiment
catalogue would require realigning the ENCODE FASTQs against CHM13.

### 11.2 Method

Coordinates are lifted through the UCSC chain file
`https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz`
(7,706 chains, 823,517 ungapped blocks, 1,273 of them minus-strand), using the
chain walker in `scripts/liftover_chm13.py` rather than the UCSC `liftOver` binary, so
that the **alignment strand is retained**. For each site, `start`, `summit` and
`end - 1` are lifted independently; where several chains cover a position the
highest-scoring one wins. `chm13_start`/`chm13_end` are the min/max of the two
lifted endpoints, so intervals stay ascending on the minus strand.

### 11.3 Columns

| column | meaning |
|---|---|
| `chm13_chrom` | CHM13 chromosome, `""` if unmapped |
| `chm13_start` / `chm13_end` | lifted interval, 0-based half-open; `-1` unless status is `unique` |
| `chm13_summit` | lifted summit; `-1` if unmapped |
| `chm13_strand` | `+` same orientation in both assemblies, `-` **opposite orientation** |
| `chm13_status` | `unique` / `gapped` / `partial` / `unmapped` (below) |

`unique` — all three positions lifted through one chain and the width is
preserved to within `GAP_TOLERANCE = 10` bp. `gapped` — both ends lifted through
one chain but across a large chain gap, so the lifted width is wrong; summit is
still usable, interval is not. `partial` — the summit lifted but an endpoint did
not, or landed in a different chain. `unmapped` — the summit has no chain.

### 11.4 Result

| | sites |
|---|---|
| mapped | 336,575 / 337,104 (99.84%) |
| `unique` | 333,952 |
| `gapped` | 2,032 |
| `partial` | 591 |
| `unmapped` | 529 |
| minus strand (orientation differs between assemblies) | 1,500 |
| lifted to a different chromosome | 205 |

### 11.5 Verification

Checks run on the output: every lifted chromosome exists in `hs1.chrom.sizes`;
every `unique` interval lies within its CHM13 chromosome bounds; every summit
lies inside its own lifted interval; `start < end` throughout; unmapped rows are
blanked. Width is preserved exactly for 95.2% of `unique` sites, within 10 bp for
99.4%.

The decisive check is **monotonicity within a chain**: for a correct lifter,
sites ordered by hg38 position must lift in strictly ascending order on a plus
chain and descending on a minus chain. Across the 114 chains carrying two or
more sites, **0 violated this**. Order is *not* monotonic when read across a
whole chromosome, because sites switch between chains — which is signal, not
error (§11.6).

### 11.6 Reference orientation flips — read this before intersecting inversions

1,500 sites lie in segments whose orientation differs between GRCh38 and CHM13.
The contiguous blocks of at least 20 sites are:

| hg38 region | Mb | sites | interpretation |
|---|---|---|---|
| chr8:8,059,041–12,488,451 | 4.43 | 610 | **8p23.1** defensin/REPD–REPP inversion |
| chr1:144,412,607–149,076,306 | 4.66 | 221 | **1q21.1** |
| chr16:21,584,025–22,437,291 | 0.85 | 131 | **16p12.2** |
| chr3:195,693,952–195,954,668 | 0.26 | 63 | 3q29 |
| chr16:32,897,720–35,909,432 | 3.01 | 60 | pericentromeric 16p11.2/16q |
| chr21:8,596,010–10,542,796 → chr14 | 1.95 | 48 | acrocentric short arm / rDNA |
| chr9:60,758,572–65,678,984 | 4.92 | 39 | pericentromeric chr9 |
| chr22:10,573,243–12,210,574 | 1.64 | 26 | acrocentric short arm |
| chr20:29,044,194–30,922,488 | 1.88 | 24 | pericentromeric |
| chr16:28,480,491–28,623,686 | 0.14 | 41 | 16p11.2 |
| chr9:40,903,945–41,508,149 | 0.60 | 21 | pericentromeric chr9 |

The first three are well-known common polymorphic inversions: the two references
simply carry opposite alleles. **This is a confound for any inversion
intersection.** "Inverted" is defined relative to a reference, so a callset built
on CHM13 and one built on GRCh38 will assign opposite orientations to the same
haplotype in these regions. Decide on one reference orientation per locus and
polarise every callset to it before merging.

The pericentromeric and acrocentric blocks are a different matter: there the
liftover itself is unreliable, and the 205 cross-chromosome lifts (mostly
subtelomeric chr1→chr16/chr17 and chr21/chr22→chr14) reflect duplication
ambiguity rather than real rearrangement. Treat cross-chromosome lifts as
suspect regardless of their status flag.
