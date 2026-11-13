"""Shared helpers for the ENCODE CTCF site pipeline.

Imported by the scripts in this folder; not meant to be run directly. See
README.md for the specification these functions implement.
"""

from __future__ import annotations

import gzip
import io
import os
import time

import numpy as np
import pandas as pd
import requests

# --------------------------------------------------------------------------- #
# ENCODE
# --------------------------------------------------------------------------- #
ENCODE = "https://www.encodeproject.org"

METADATA_URL = (
    ENCODE + "/metadata/?type=Experiment&assay_title=TF+ChIP-seq"
    "&target.label=CTCF&assembly=GRCh38&status=released"
    "&files.file_type=bed+narrowPeak&files.assembly=GRCh38"
)

# ENCODE4 uniform-pipeline output preferred, then legacy ENCODE3 flavours.
OUTPUT_RANK = {
    "IDR thresholded peaks": 0,
    "optimal IDR thresholded peaks": 1,
    "conservative IDR thresholded peaks": 2,
    "pseudoreplicated IDR thresholded peaks": 3,
}
ANALYSIS_RANK = {
    "ENCODE4 v1.8.0 GRCh38": 0,
    "ENCODE4 v1.7.0 GRCh38": 1,
    "ENCODE4 v1.6.1 GRCh38": 2,
    "ENCODE4 v1.5.1 GRCh38": 3,
    "ENCODE4 v1.5.0 GRCh38": 4,
    "ENCODE4 v1.4.0 GRCh38": 5,
    "ENCODE3 GRCh38": 6,
    "Lab custom GRCh38": 7,
}
# Audits meaning the peak set itself is untrustworthy.
FATAL_AUDITS = ("extremely low read depth", "control extremely low read depth")

MANIFEST_COLUMNS = [
    "experiment", "file_accession", "biosample", "biosample_type", "output_type",
    "analysis", "lab", "date_released", "size_bytes", "md5sum", "url",
    "audit_error", "audit_not_compliant",
]

# --------------------------------------------------------------------------- #
# Assembly
# --------------------------------------------------------------------------- #
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CHROM_CODE = {c: i for i, c in enumerate(CHROMS)}

SLOP = 75          # half-width of the summit window; summits <= 2*SLOP apart merge
GAP_TOLERANCE = 10  # max |lifted width - hg38 width| still called a clean lift

CHAIN_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz"
)
HG38_CHROM_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{chrom}.fa.gz"
CHR21_URL = HG38_CHROM_URL.format(chrom="chr21")
JASPAR_URL = "https://jaspar.elixir.no/api/v1/matrix/MA0139.1/?format=json"
HS1_SIZES_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.chrom.sizes"
# JASPAR's own genome-wide MA0139.1 predictions for hg38 (JASPAR 2022 UCSC
# tracks; later releases carry only the successor matrix MA0139.2)
JASPAR_TRACK_URL = (
    "https://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2022/hg38/MA0139.1.tsv.gz"
)
JASPAR_TRACK_COLUMNS = ["chrom", "start", "end", "name", "rel_score_1000",
                        "pvalue_score", "strand"]

# motif scan (README §9): a window is a hit at or above this JASPAR relative
# score; the best hit is searched within +/-MOTIF_FLANK of the summit
MOTIF_REL_CUTOFF = 0.80
MOTIF_FLANK = 100


# --------------------------------------------------------------------------- #
def fetch(url: str, retries: int = 4, timeout: int = 600) -> bytes:
    """GET with exponential backoff."""
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=timeout)
            r.raise_for_status()
            return r.content
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def select_files(meta: pd.DataFrame, include_low_quality: bool = False) -> pd.DataFrame:
    """One peak file per experiment: best output type from the newest analysis."""
    f = meta[
        meta["Output type"].isin(OUTPUT_RANK)
        & (meta["File Status"] == "released")
        & (meta["File analysis status"] == "released")
        & (meta["File assembly"] == "GRCh38")
        & (meta["Biosample organism"] == "Homo sapiens")
    ].copy()

    f["output_rank"] = f["Output type"].map(OUTPUT_RANK)
    f["analysis_rank"] = f["File analysis title"].map(ANALYSIS_RANK).fillna(99)
    f = f.sort_values(
        ["Experiment accession", "analysis_rank", "output_rank", "Size"],
        ascending=[True, True, True, False],
    )
    sel = f.groupby("Experiment accession", as_index=False).first()

    audit = sel["Audit ERROR"].fillna("")
    sel["fatal_audit"] = audit.str.contains("|".join(FATAL_AUDITS), case=False)
    sel["audit_error"] = audit
    sel["audit_not_compliant"] = sel["Audit NOT_COMPLIANT"].fillna("")
    n_bad = int(sel["fatal_audit"].sum())
    if not include_low_quality:
        sel = sel[~sel["fatal_audit"]]
    sel.attrs["n_dropped"] = n_bad

    sel = sel.rename(columns={
        "File accession": "file_accession", "Experiment accession": "experiment",
        "Output type": "output_type", "Biosample term name": "biosample",
        "Biosample type": "biosample_type", "File analysis title": "analysis",
        "File download URL": "url", "Size": "size_bytes", "Lab": "lab",
        "Experiment date released": "date_released",
    })
    out = sel[MANIFEST_COLUMNS].reset_index(drop=True)
    out.attrs["n_dropped"] = n_bad
    return out


def biosample_codes(manifest: pd.DataFrame) -> dict[str, int]:
    """Stable biosample -> index map, identical for every chunk."""
    return {b: i for i, b in enumerate(sorted(manifest["biosample"].unique()))}


def parse_peaks(raw: bytes):
    """narrowPeak bytes -> (chrom_code, summit, signalValue, -log10 qvalue)."""
    df = pd.read_csv(
        io.BytesIO(gzip.decompress(raw)),
        sep="\t", header=None, usecols=[0, 1, 2, 6, 8, 9],
        names=["chrom", "start", "end", "signal", "negl10q", "offset"],
        dtype={"chrom": "string", "start": np.int64, "end": np.int64,
               "signal": np.float32, "negl10q": np.float32, "offset": np.int64},
    )
    code = df["chrom"].map(CHROM_CODE).astype("float64").to_numpy()
    keep = ~np.isnan(code)
    df = df[keep]
    code = code[keep].astype(np.int8)

    start = df["start"].to_numpy()
    offset = df["offset"].to_numpy()
    mid = (start + df["end"].to_numpy()) // 2
    summit = np.where(offset >= 0, start + offset, mid).astype(np.int32)
    return code, summit, df["signal"].to_numpy(np.float32), df["negl10q"].to_numpy(np.float32)


def load_chunks(paths: list[str]) -> dict:
    """Concatenate per-chunk .npz peak arrays into one pooled set."""
    parts = {k: [] for k in ("chrom", "summit", "signal", "negq", "exp", "bio")}
    for p in sorted(paths):
        z = np.load(p)
        for k in parts:
            parts[k].append(z[k])
    return {k: np.concatenate(v) for k, v in parts.items()}


def build_consensus(P: dict, n_exp_total: int, n_bio_total: int,
                    slop: int = SLOP) -> pd.DataFrame:
    """Cluster summits (single linkage, gap <= 2*slop) into consensus sites."""
    order = np.lexsort((P["summit"], P["chrom"]))
    chrom = P["chrom"][order]
    summit = P["summit"][order].astype(np.int64)
    signal = P["signal"][order]
    negq = P["negq"][order]
    exp = P["exp"][order].astype(np.int64)
    bio = P["bio"][order].astype(np.int64)

    brk = np.ones(summit.size, dtype=bool)
    brk[1:] = (chrom[1:] != chrom[:-1]) | (np.diff(summit) > 2 * slop)
    starts = np.flatnonzero(brk)
    cid = np.cumsum(brk) - 1
    n_clusters = starts.size
    ends = np.append(starts[1:], summit.size)

    max_signal = np.maximum.reduceat(signal, starts)
    sum_signal = np.add.reduceat(signal.astype(np.float64), starts)
    max_negq = np.maximum.reduceat(negq, starts)
    n_peaks = ends - starts

    # representative summit = summit of the strongest contributing peak
    is_max = signal == max_signal[cid]
    idx = np.where(is_max, np.arange(summit.size), summit.size)
    rep_summit = summit[np.minimum.reduceat(idx, starts)]

    n_exp = np.bincount(np.unique(cid * n_exp_total + exp) // n_exp_total,
                        minlength=n_clusters)
    n_bio = np.bincount(np.unique(cid * n_bio_total + bio) // n_bio_total,
                        minlength=n_clusters)

    lo = summit[starts] - slop
    hi = summit[ends - 1] + slop + 1

    sites = pd.DataFrame({
        "chrom": pd.Categorical.from_codes(chrom[starts], categories=CHROMS),
        "start": np.maximum(lo, 0).astype(np.int32),
        "end": hi.astype(np.int32),
        "summit": rep_summit.astype(np.int32),
        "n_experiments": n_exp.astype(np.int16),
        "n_biosamples": n_bio.astype(np.int16),
        "frac_experiments": (n_exp / n_exp_total).astype(np.float32),
        "max_signal": max_signal,
        "mean_signal": (sum_signal / n_peaks).astype(np.float32),
        "max_neglog10q": max_negq,
        "n_peaks": n_peaks.astype(np.int32),
    })
    sites["width"] = (sites["end"] - sites["start"]).astype(np.int32)
    return sites.sort_values(["chrom", "start"], ignore_index=True)


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #
DTYPES = {
    "chrom": "category", "start": np.int32, "end": np.int32, "summit": np.int32,
    "n_experiments": np.int16, "n_biosamples": np.int16,
    "frac_experiments": np.float32, "max_signal": np.float32,
    "mean_signal": np.float32, "max_neglog10q": np.float32,
    "n_peaks": np.int32, "width": np.int32,
    # motif orientation (README §9.2)
    "motif_rel": np.float32, "motif_strand": "category", "motif_start": np.int32,
    "n_motifs_plus": np.int16, "n_motifs_minus": np.int16,
    # present only in the dual-coordinate table
    "chm13_chrom": "string", "chm13_start": np.int64, "chm13_end": np.int64,
    "chm13_summit": np.int64, "chm13_strand": "string",
    "motif_strand_chm13": "category", "chm13_status": "category",
}


def load_ctcf_sites(path: str, min_experiments: int = 1) -> pd.DataFrame:
    """Load a CTCF site table, filtering on cross-experiment recurrence.

    Accepts either a parquet dataset written by `pd_lfs.write_parquet` (a
    directory, or an https URL to one) or a plain `.tsv.gz`. Parquet carries its
    own dtypes; the TSV path applies DTYPES explicitly.
    """
    if str(path).startswith(("http://", "https://")) or os.path.isdir(path):
        from pd_lfs import read_parquet
        df = read_parquet(path)
    else:
        header = pd.read_csv(path, sep="\t", nrows=0).columns
        dt = {k: v for k, v in DTYPES.items() if k in header}
        df = pd.read_csv(path, sep="\t", dtype=dt)
    if min_experiments > 1:
        df = df[df["n_experiments"] >= min_experiments].reset_index(drop=True)
    return df


# --------------------------------------------------------------------------- #
# hg38 -> CHM13 liftover
# --------------------------------------------------------------------------- #
class Chain:
    __slots__ = ("score", "t_name", "q_name", "q_size", "q_rev",
                 "t_start", "t_end", "q_start", "t_lo", "t_hi")


def load_chains(path: str) -> list:
    """Parse a UCSC .over.chain file into per-chain ungapped block arrays."""
    chains: list = []
    cur = None
    tb: list = []
    qb: list = []
    sz: list = []
    tpos = qpos = 0

    def close():
        if cur is None or not sz:
            return
        cur.t_start = np.asarray(tb, np.int64)
        cur.q_start = np.asarray(qb, np.int64)
        cur.t_end = cur.t_start + np.asarray(sz, np.int64)
        cur.t_lo, cur.t_hi = int(cur.t_start[0]), int(cur.t_end[-1])
        chains.append(cur)

    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("chain"):
                close()
                f = line.split()
                cur = Chain()
                cur.score = int(f[1])
                cur.t_name, cur.q_name = f[2], f[7]
                cur.q_size = int(f[8])
                cur.q_rev = f[9] == "-"
                tpos, qpos = int(f[5]), int(f[10])
                tb, qb, sz = [], [], []
            else:
                f = line.split()
                size = int(f[0])
                tb.append(tpos); qb.append(qpos); sz.append(size)
                if len(f) == 3:
                    tpos += size + int(f[1])
                    qpos += size + int(f[2])
    close()
    return chains


def lift(chains: list, chrom: np.ndarray, pos: np.ndarray):
    """Lift (chrom, pos) through the chains; the best-scoring chain wins.

    Returns (q_chrom, q_pos, strand, chain_id); unmapped entries are ""/-1/0/-1.
    """
    n = pos.size
    q_chrom = np.full(n, "", dtype=object)
    q_pos = np.full(n, -1, np.int64)
    strand = np.zeros(n, np.int8)
    chain_id = np.full(n, -1, np.int32)

    by_chrom: dict = {}
    for c in np.unique(chrom):
        idx = np.flatnonzero(chrom == c)
        o = np.argsort(pos[idx], kind="stable")
        by_chrom[c] = (idx[o], pos[idx][o])

    # ascending score, so a higher-scoring chain overwrites a lower-scoring one
    for ch in sorted(range(len(chains)), key=lambda i: chains[i].score):
        C = chains[ch]
        hit = by_chrom.get(C.t_name)
        if hit is None:
            continue
        idx_sorted, p_sorted = hit
        lo = np.searchsorted(p_sorted, C.t_lo, "left")
        hi = np.searchsorted(p_sorted, C.t_hi, "right")
        if lo >= hi:
            continue
        p = p_sorted[lo:hi]
        j = np.searchsorted(C.t_start, p, "right") - 1
        ok = (j >= 0) & (p < C.t_end[np.maximum(j, 0)])
        if not ok.any():
            continue
        j = j[ok]
        target = idx_sorted[lo:hi][ok]
        raw = C.q_start[j] + (p[ok] - C.t_start[j])
        q_chrom[target] = C.q_name
        q_pos[target] = (C.q_size - 1 - raw) if C.q_rev else raw
        strand[target] = -1 if C.q_rev else 1
        chain_id[target] = ch
    return q_chrom, q_pos, strand, chain_id


# --------------------------------------------------------------------------- #
# CTCF motif scanning (orientation annotation and validation)
# --------------------------------------------------------------------------- #
def load_pwm(pfm_json: dict):
    """JASPAR PFM -> (log-odds, reverse-complement log-odds, min score, max score)."""
    M = np.array([pfm_json[b] for b in "ACGT"], float)
    pseudo = 0.01 * M.sum(0).mean()
    freq = (M + 0.25 * pseudo) / (M.sum(0) + pseudo)
    lo = np.log2(freq / 0.25)
    return lo, lo[::-1, ::-1], lo.min(0).sum(), lo.max(0).sum()


def encode_sequence(seq: str) -> np.ndarray:
    """A,C,G,T -> 0..3; everything else -> -1."""
    code = np.full(256, -1, np.int8)
    for i, b in enumerate("ACGT"):
        code[ord(b)] = i
    return code[np.frombuffer(seq.upper().encode(), dtype=np.uint8)]


def best_relscore(S, centers, lo, rc, mn, mx, flank: int = 100) -> np.ndarray:
    """Max relative PWM score on either strand within +/-flank of each center."""
    W = lo.shape[1]
    win = S[centers[:, None] + np.arange(-flank, flank + 1)[None, :]]
    bad = win < 0
    win = np.where(bad, 0, win)
    npos = win.shape[1] - W + 1
    fwd = np.zeros((len(centers), npos))
    rev = np.zeros((len(centers), npos))
    nmask = np.zeros((len(centers), npos), bool)
    for j in range(W):
        col = win[:, j:j + npos]
        fwd += lo[col, j]
        rev += rc[col, j]
        nmask |= bad[:, j:j + npos]
    sc = np.maximum(fwd, rev)
    sc[nmask] = -np.inf
    return (sc.max(1) - mn) / (mx - mn)


def scan_pwm(S: np.ndarray, lo, rc, mn, mx, block: int = 1 << 24):
    """Relative PWM score of the W-mer starting at every position, both strands.

    Returns `(fwd, rev)` float32 arrays of length `len(S) - W + 1`, holding the
    JASPAR relative score `(score - min) / (max - min)` of the window starting
    at each position on the plus strand and on the minus strand respectively.
    Windows containing a non-ACGT base score -1 (real scores lie in [0, 1]).
    The chromosome is processed in blocks so that the temporaries stay small.
    """
    W = lo.shape[1]
    n = S.size - W + 1
    lo32, rc32 = lo.astype(np.float32), rc.astype(np.float32)
    scale = np.float32(mx - mn)
    fwd = np.empty(n, np.float32)
    rev = np.empty(n, np.float32)
    for b0 in range(0, n, block):
        b1 = min(b0 + block, n)
        m = b1 - b0
        seg = S[b0:b1 + W - 1]
        bad = seg < 0
        seg = np.where(bad, 0, seg).astype(np.intp)
        f = np.zeros(m, np.float32)
        r = np.zeros(m, np.float32)
        nmask = np.zeros(m, bool)
        for j in range(W):
            col = seg[j:j + m]
            f += lo32[col, j]
            r += rc32[col, j]
            nmask |= bad[j:j + m]
        f = (f - np.float32(mn)) / scale
        r = (r - np.float32(mn)) / scale
        f[nmask] = -1
        r[nmask] = -1
        fwd[b0:b1] = f
        rev[b0:b1] = r
    return fwd, rev


def annotate_motifs(sites: pd.DataFrame, fwd: np.ndarray, rev: np.ndarray, width: int,
                    flank: int = MOTIF_FLANK, cutoff: float = MOTIF_REL_CUTOFF) -> pd.DataFrame:
    """Per-site motif orientation from the strand score arrays of one chromosome.

    `sites` holds the rows of one chromosome (`start`, `end`, `summit`); `fwd`
    and `rev` come from `scan_pwm` on that chromosome. The best window on either
    strand within +/-`flank` of the summit (the same window `best_relscore` uses)
    gives `motif_rel`; if it reaches `cutoff` the site is oriented by that hit
    (`motif_strand`, `motif_start`), otherwise `motif_strand` is `.` and
    `motif_start` is -1. `n_motifs_plus`/`n_motifs_minus` count every window at
    or above `cutoff` on each strand whose start lies inside the whole site
    interval, so sites carrying motifs in both orientations can be recognised.
    """
    n = fwd.size
    summit = sites["summit"].to_numpy(np.int64)
    offs = np.arange(-flank, flank - width + 2)               # window starts
    idx = summit[:, None] + offs[None, :]
    valid = (idx >= 0) & (idx < n)
    idx = np.clip(idx, 0, n - 1)
    f = np.where(valid, fwd[idx], np.float32(-1))
    r = np.where(valid, rev[idx], np.float32(-1))
    best_f, best_r = f.max(1), r.max(1)
    is_plus = best_f >= best_r
    rel = np.where(is_plus, best_f, best_r)
    start = summit + offs[np.where(is_plus, f.argmax(1), r.argmax(1))]
    hit = rel >= cutoff

    # hits are sparse, so count them per site from their sorted positions
    hp = np.flatnonzero(fwd >= cutoff)
    hm = np.flatnonzero(rev >= cutoff)
    a = sites["start"].to_numpy(np.int64)
    b = np.maximum(a, sites["end"].to_numpy(np.int64) - width + 1)   # exclusive

    return pd.DataFrame({
        "motif_rel": rel.astype(np.float32),
        "motif_strand": np.where(hit, np.where(is_plus, "+", "-"), "."),
        "motif_start": np.where(hit, start, -1).astype(np.int32),
        "n_motifs_plus": (np.searchsorted(hp, b) - np.searchsorted(hp, a)).astype(np.int16),
        "n_motifs_minus": (np.searchsorted(hm, b) - np.searchsorted(hm, a)).astype(np.int16),
    }, index=sites.index)


def read_fasta_gz(path: str) -> str:
    """Sequence of a single-record gzipped FASTA, as one string."""
    with gzip.open(path, "rt") as fh:
        return "".join(l.strip() for l in fh if not l.startswith(">"))
