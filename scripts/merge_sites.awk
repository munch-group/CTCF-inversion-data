# One-pass consensus builder -- a dependency-free reimplementation of the
# clustering and annotation steps (METHODS.md §4-5, §7).
#
# Input : summit BED sorted with `LC_ALL=C sort -k1,1 -k2,2n`, tab-separated:
#         chrom, summit, summit+1, experiment, biosample, signalValue, -log10(q)
# Output: chrom, start, end, summit, n_experiments, n_biosamples, n_peaks,
#         max_signal, mean_signal, max_neglog10q
# Run as: LC_ALL=C awk -f merge_sites.awk summits.sorted.bed
BEGIN { FS="\t"; OFS="\t"; if (SLOP=="") SLOP=75; first=1 }
function emit() {
    st = cstart - SLOP; if (st < 0) st = 0
    printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\n", \
        pchrom, st, clast + SLOP + 1, bestsummit, nexp, nbio, npk, maxsig, sumsig/npk, maxq
}
{
    chrom=$1; s=$2+0; ex=$4; bio=$5; sig=$6+0; q=$7+0
    if (first || chrom != pchrom || s - clast > 2*SLOP) {
        if (!first) emit()
        cstart=s; nexp=0; nbio=0; npk=0; sumsig=0; maxsig=-1; maxq=-1
        split("", E); split("", B); first=0
    }
    npk++; sumsig += sig
    if (sig > maxsig) { maxsig=sig; bestsummit=s }
    if (q > maxq) maxq=q
    if (!(ex in E)) { E[ex]=1; nexp++ }
    if (!(bio in B)) { B[bio]=1; nbio++ }
    pchrom=chrom; clast=s
}
END { if (!first) emit() }
