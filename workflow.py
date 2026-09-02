# %% [markdown]
# ---
# title: GWF workflow
# execute:
#   eval: false
# ---

# %% [markdown]
r"""
Retrieval, compilation and validation of the ENCODE CTCF binding-site catalogue
in hg38, with T2T-CHM13v2.0 coordinates added by liftover.

Every released human CTCF TF ChIP-seq experiment on GRCh38 is pulled from the
ENCODE portal, one IDR-thresholded peak file per experiment is selected, the
peak summits are pooled and clustered into consensus sites, and the result is
checked against independent sequence evidence. See `METHODS.md` for the full
specification — the section numbers referenced in the templates below point
into it.

```{mermaid}
flowchart TD
    A[fetch_encode_metadata<br/><i>steps/encode</i>] --> B[select_encode_files<br/><b>results/ctcf_encode_files.tsv</b>]
    B --> C[fetch_peaks_0 .. fetch_peaks_15<br/><i>steps/peaks</i>]
    C --> D[build_consensus<br/><i>steps/consensus</i>]
    R1[fetch_reference chain] --> E
    D --> E[liftover_chm13<br/><i>steps/consensus</i>]
    D --> P1[export_sites_parquet<br/><b>results/ctcf_sites_hg38/</b>]
    E --> P2[export_sites_chm13_parquet<br/><b>results/ctcf_sites_hg38_chm13/</b>]
    D --> V1[validate_motif]
    R2[fetch_reference chr21] --> V1
    R3[fetch_reference jaspar] --> V1
    E --> V2[validate_liftover]
    R4[fetch_reference hs1_sizes] --> V2
    R1 --> V2
    C --> V3[validate_consensus<br/><i>via merge_sites.awk</i>]
    D --> V3
    P1 --> N[notebooks]
    P2 --> N
    V1 --> N[notebooks]
    V2 --> N
    V3 --> N
```

Files under `steps/` are intermediates and are git-ignored; everything worth
keeping is written to `results/`. The two large site tables are published there
as parquet datasets rather than TSV, sharded into part files below GitHub's
50 MB limit and readable over plain HTTPS with `pd_lfs.read_parquet`.
"""

# %% [markdown]
"""
## Imports and utility functions
"""

# %%
import glob
import os
import re
from pathlib import Path

from gwf import AnonymousTarget, Workflow
from gwf.workflow import collect

# directories
STEPS = 'steps'                     # intermediate / temporary files (git-ignored)
RESULTS = 'results'                 # files worth keeping
TMP = f'{STEPS}/tmp'                # write-then-move scratch, same filesystem as results

# number of parallel download/parse tasks the peak files are split over. Fixed
# here rather than derived from the manifest, because gwf must know the shape of
# the graph when this file is evaluated -- before the manifest exists.
CHUNKS = 16


# %%
def modify_path(path, **kwargs):
    """
    Utility function for modifying file paths substituting
    the directory (dir), base name (base), or file suffix (suffix).
    """
    for key in ['dir', 'base', 'suffix']:
        kwargs.setdefault(key, None)
    assert len(kwargs) == 3

    par, name = os.path.split(path)
    name_no_suffix, suf = os.path.splitext(name)
    if type(kwargs['suffix']) is str:
        suf = kwargs['suffix']
    if kwargs['dir'] is not None:
        par = kwargs['dir']
    if kwargs['base'] is not None:
        name_no_suffix = kwargs['base']

    new_path = os.path.join(par, name_no_suffix + suf)
    if type(kwargs['suffix']) is tuple:
        assert len(kwargs['suffix']) == 2
        new_path, nsubs = re.subn(r'{}$'.format(kwargs['suffix'][0]), kwargs['suffix'][1], new_path)
        assert nsubs == 1, nsubs
    return new_path


def tmp_path(path):
    """Scratch path under steps/ mirroring the basename of an output file."""
    return os.path.join(TMP, os.path.basename(path))


# %% [markdown]
"""
## Template functions

Each task runs through `pixi run`, so it executes in the project environment
regardless of the environment the worker pool was started in. Plain POSIX
utilities (`mkdir`, `mv`, `sort`, `awk`) are invoked directly -- they come from
the system, not from the pixi environment.

Every template writes to a scratch file under `steps/tmp` and moves it into
place only if the command succeeded, so a crash can never leave a partial
output that GWF would mistake for a finished one. The validation templates rely
on the same idiom for a second purpose: their scripts exit non-zero when a check
fails, which stops the `&&` chain, so a failed validation produces no report and
the target stays incomplete rather than silently passing.
"""

# %% [markdown]
"""
### Retrieval
"""


# %%
def fetch_encode_metadata():
    """
    Downloads the ENCODE batch-metadata table for CTCF narrowPeak files (METHODS §2.1).
    """
    output_dir = f'{STEPS}/encode'
    metadata_path = f'{output_dir}/encode_ctcf_metadata.tsv'

    inputs = []
    outputs = {'metadata_path': metadata_path}
    options = {'memory': '2g', 'walltime': '01:00:00'}

    tmp = tmp_path(metadata_path)
    spec = f"""
    mkdir -p {output_dir} {TMP}
    pixi run python scripts/fetch_encode_metadata.py --out {tmp} &&
        mv {tmp} {metadata_path}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def select_encode_files(metadata_path):
    """
    Picks one peak file per experiment and writes the provenance manifest (METHODS §2.2-2.3).
    """
    manifest_path = f'{RESULTS}/ctcf_encode_files.tsv'

    inputs = [metadata_path]
    outputs = {'manifest_path': manifest_path}
    options = {'memory': '4g', 'walltime': '00:20:00'}

    tmp = tmp_path(manifest_path)
    spec = f"""
    mkdir -p {RESULTS} {TMP}
    pixi run python scripts/select_encode_files.py --metadata {metadata_path} --out {tmp} &&
        mv {tmp} {manifest_path}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def fetch_reference(what):
    """
    Downloads one external reference file (liftover chains, chr21, JASPAR motif, CHM13 sizes).
    """
    output_dir = f'{STEPS}/refs'
    suffix = {'chain': 'hg38ToHs1.over.chain.gz', 'chr21': 'hg38_chr21.fa.gz',
              'jaspar': 'MA0139.1.json', 'hs1_sizes': 'hs1.chrom.sizes'}[what]
    ref_path = f'{output_dir}/{suffix}'

    inputs = []
    outputs = {'ref_path': ref_path}
    options = {'memory': '2g', 'walltime': '02:00:00'}

    tmp = tmp_path(ref_path)
    spec = f"""
    mkdir -p {output_dir} {TMP}
    pixi run python scripts/fetch_reference.py --what {what} --out {tmp} &&
        mv {tmp} {ref_path}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def fetch_peak_chunk(manifest_path, chunk, n_chunks):
    """
    Downloads and parses one chunk of the ENCODE peak files into summit arrays (METHODS §3).
    """
    output_dir = f'{STEPS}/peaks'
    chunk_path = f'{output_dir}/chunk_{chunk:02d}.npz'

    inputs = [manifest_path]
    outputs = {'chunk_path': chunk_path}
    options = {'cores': 4, 'memory': '8g', 'walltime': '04:00:00'}

    tmp = tmp_path(chunk_path)
    spec = f"""
    mkdir -p {output_dir} {TMP}
    pixi run python scripts/fetch_peak_chunk.py \\
        --manifest {manifest_path} --chunk {chunk} --chunks {n_chunks} \\
        --workers 4 --out {tmp} &&
        mv {tmp} {chunk_path}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# %% [markdown]
"""
### Compilation
"""


# %%
def build_consensus(manifest_path, chunk_paths):
    """
    Clusters the pooled summits into the consensus CTCF site table (METHODS §4-5).
    """
    output_dir = f'{STEPS}/consensus'
    sites_path = f'{output_dir}/ctcf_sites_hg38.tsv.gz'

    inputs = [manifest_path] + list(chunk_paths)
    outputs = {'sites_path': sites_path}
    options = {'memory': '32g', 'walltime': '02:00:00'}

    tmp = tmp_path(sites_path)
    spec = f"""
    mkdir -p {output_dir} {TMP}
    pixi run python scripts/build_consensus.py \\
        --manifest {manifest_path} --chunks {' '.join(chunk_paths)} --out {tmp} &&
        mv {tmp} {sites_path}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def liftover_chm13(sites_path, chain_path):
    """
    Appends T2T-CHM13v2.0 coordinates, keeping the alignment strand (METHODS §11).
    """
    output_dir = f'{STEPS}/consensus'
    dual_path = f'{output_dir}/ctcf_sites_hg38_chm13.tsv.gz'

    inputs = [sites_path, chain_path]
    outputs = {'dual_path': dual_path}
    options = {'memory': '16g', 'walltime': '01:00:00'}

    tmp = tmp_path(dual_path)
    spec = f"""
    mkdir -p {output_dir} {TMP}
    pixi run python scripts/liftover_chm13.py \\
        --sites {sites_path} --chain {chain_path} --out {tmp} &&
        mv {tmp} {dual_path}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def export_parquet(tsv_path, name):
    """
    Converts a site table into a parquet dataset of GitHub-sized part files.

    `pd_lfs.write_parquet` shards the frame into `part-*.parquet` files each
    below 40 MB -- comfortably under GitHub's 50 MB warning -- and writes a
    `_manifest.json` index that lets `pd_lfs.read_parquet` read the dataset over
    plain HTTPS. The manifest is the tracked output: it exists only once the
    whole dataset has been written.
    """
    dataset_dir = f'{RESULTS}/{name}'
    manifest = f'{dataset_dir}/_manifest.json'

    inputs = [tsv_path]
    outputs = {'manifest_path': manifest}
    options = {'memory': '16g', 'walltime': '00:30:00'}

    # a dataset is a directory, so the write-then-move idiom moves the whole
    # directory; the old one is only removed once the new one is complete
    tmp_dir = f'{TMP}/{name}'
    spec = f"""
    mkdir -p {RESULTS} {TMP}
    rm -rf {tmp_dir} &&
    pixi run python scripts/export_parquet.py --sites {tsv_path} --out {tmp_dir} &&
        rm -rf {dataset_dir} && mv {tmp_dir} {dataset_dir}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# %% [markdown]
"""
### Validation

Three independent checks. `validate_motif` is the strongest: it uses sequence
evidence that never entered the compilation, so it confirms both the assembly
and the biological identity of the sites. `validate_consensus` re-derives one
chromosome with a dependency-free awk implementation of the same specification
and diffs it against the Python result. `validate_liftover` checks the derived
CHM13 coordinates and reports the loci where the two references disagree on
orientation.
"""


# %%
def validate_motif(sites_path, chr21_path, pfm_path):
    """
    Scans the JASPAR CTCF motif around chr21 summits against shuffled controls (METHODS §9).
    """
    output_dir = f'{RESULTS}/validation'
    report_path = f'{output_dir}/motif_enrichment.txt'

    inputs = [sites_path, chr21_path, pfm_path]
    outputs = {'report_path': report_path}
    options = {'memory': '8g', 'walltime': '00:30:00'}

    tmp = tmp_path(report_path)
    spec = f"""
    mkdir -p {output_dir} {TMP}
    pixi run python scripts/validate_motif.py \\
        --sites {sites_path} --chr21 {chr21_path} --pfm {pfm_path} --out {tmp} &&
        mv {tmp} {report_path}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def validate_liftover(dual_path, chain_path, sizes_path):
    """
    Checks the CHM13 coordinates and reports reference orientation flips (METHODS §11.5-11.6).
    """
    output_dir = f'{RESULTS}/validation'
    report_path = f'{output_dir}/liftover_checks.txt'

    inputs = [dual_path, chain_path, sizes_path]
    outputs = {'report_path': report_path}
    options = {'memory': '16g', 'walltime': '00:30:00'}

    tmp = tmp_path(report_path)
    spec = f"""
    mkdir -p {output_dir} {TMP}
    pixi run python scripts/validate_liftover.py \\
        --sites {dual_path} --chain {chain_path} --chrom-sizes {sizes_path} --out {tmp} &&
        mv {tmp} {report_path}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def validate_consensus(manifest_path, chunk_paths, sites_path, chrom='chr21'):
    """
    Re-derives one chromosome with merge_sites.awk and diffs it against the table (METHODS §7).
    """
    work_dir = f'{STEPS}/validation'
    output_dir = f'{RESULTS}/validation'
    report_path = f'{output_dir}/consensus_awk_check.txt'

    bed = f'{work_dir}/{chrom}_summits.bed'
    sorted_bed = f'{work_dir}/{chrom}_summits.sorted.bed'
    awk_sites = f'{work_dir}/{chrom}_awk_sites.tsv'

    inputs = [manifest_path, sites_path] + list(chunk_paths)
    outputs = {'report_path': report_path}
    options = {'memory': '16g', 'walltime': '00:30:00'}

    tmp = tmp_path(report_path)
    # LC_ALL=C throughout: a comma-decimal locale corrupts awk's numeric output
    # and locale collation changes sort order.
    spec = f"""
    mkdir -p {work_dir} {output_dir} {TMP}
    pixi run python scripts/dump_summits.py \\
        --manifest {manifest_path} --chunks {' '.join(chunk_paths)} \\
        --chrom {chrom} --out {bed} &&
    LC_ALL=C sort -k1,1 -k2,2n {bed} > {sorted_bed} &&
    LC_ALL=C awk -f scripts/merge_sites.awk {sorted_bed} > {awk_sites} &&
    pixi run python scripts/compare_consensus.py \\
        --awk {awk_sites} --sites {sites_path} --chrom {chrom} --out {tmp} &&
        mv {tmp} {report_path}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# %% [markdown]
"""
### Documentation
"""


# %%
def run_notebook(path, dependencies, memory='8g', walltime='00:30:00', cores=1):
    """
    Executes a notebook inplace and saves the output.
    """
    sentinel = modify_path(path, base=f'.{str(Path(path).name)}', suffix='.sentinel')

    inputs = [path] + dependencies
    outputs = {'sentinel': sentinel}
    options = {'memory': memory, 'walltime': walltime, 'cores': cores}

    spec = f"""
    pixi run jupyter nbconvert --to notebook --execute --inplace {path} && touch {sentinel}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# %% [markdown]
"""
## Workflow
"""

# %%
gwf = Workflow(defaults={'account': 'CTCF-inversion-data'})

# --- retrieval ------------------------------------------------------------- #

metadata_target = gwf.target_from_template(
    'fetch_encode_metadata', fetch_encode_metadata())

manifest_target = gwf.target_from_template(
    'select_encode_files',
    select_encode_files(metadata_target.outputs['metadata_path']))

manifest_path = manifest_target.outputs['manifest_path']

# external reference files, one target each
reference_targets = gwf.map(
    fetch_reference,
    [{'what': w} for w in ['chain', 'chr21', 'jaspar', 'hs1_sizes']],
    name=lambda idx, target: f'fetch_reference_{target.spec.split("--what ")[1].split()[0]}',
)
chain_path, chr21_path, pfm_path, sizes_path = [
    o['ref_path'] for o in reference_targets.outputs]

# the peak files, split over a fixed number of parallel tasks
peak_targets = gwf.map(
    fetch_peak_chunk,
    [{'chunk': i} for i in range(CHUNKS)],
    extra=dict(manifest_path=manifest_path, n_chunks=CHUNKS),
    name='fetch_peaks',
)
chunk_paths = collect(peak_targets.outputs, ['chunk_path'])['chunk_paths']

# --- compilation ----------------------------------------------------------- #

sites_target = gwf.target_from_template(
    'build_consensus', build_consensus(manifest_path, chunk_paths))
sites_path = sites_target.outputs['sites_path']

liftover_target = gwf.target_from_template(
    'liftover_chm13', liftover_chm13(sites_path, chain_path))
dual_path = liftover_target.outputs['dual_path']

# the large tables are published as parquet datasets, not as TSV -- see the
# export_parquet template for why
gwf.target_from_template(
    'export_sites_parquet', export_parquet(sites_path, 'ctcf_sites_hg38'))

gwf.target_from_template(
    'export_sites_chm13_parquet',
    export_parquet(dual_path, 'ctcf_sites_hg38_chm13'))

# --- validation ------------------------------------------------------------ #

gwf.target_from_template(
    'validate_motif', validate_motif(sites_path, chr21_path, pfm_path))

gwf.target_from_template(
    'validate_liftover', validate_liftover(dual_path, chain_path, sizes_path))

gwf.target_from_template(
    'validate_consensus', validate_consensus(manifest_path, chunk_paths, sites_path))

# --- documentation --------------------------------------------------------- #

# make notebooks depend on all output files from workflow
notebook_dependencies = []
for x in gwf.targets.values():
    outputs = x.outputs
    if type(outputs) is dict:
        for o in outputs.values():
            notebook_dependencies.append(o)
    elif type(outputs) is list:
        notebook_dependencies.extend(outputs)

#  run notebooks in sorted order nb01_, nb02_, ...
for path in sorted(glob.glob('notebooks/*.ipynb')):
    target = gwf.target_from_template(
        os.path.basename(path), run_notebook(path, notebook_dependencies))
    # make notebooks depend on all previous notebooks
    notebook_dependencies.append(target.outputs['sentinel'])

# %%
