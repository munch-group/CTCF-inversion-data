"""Shared definitions for the T2T ape all-to-all chain analysis.

The alignment is the 8-way Cactus alignment of T2T ape assemblies published at
https://cgl.gi.ucsc.edu/data/cactus/t2t-apes/8-t2t-apes-2023v2/, exported to
all-to-all chains with `cactus-hal2chains`. See README.md §12.
"""

from __future__ import annotations

CHAIN_BASE = ("https://cgl.gi.ucsc.edu/data/cactus/t2t-apes/"
              "8-t2t-apes-2023v2/chains")

# accession -> readable species name (from the alignment README)
SPECIES = {
    'GCA_028858775.2': 'chimpanzee',           # Pan troglodytes,      mPanTro3
    'GCA_029289425.2': 'bonobo',               # Pan paniscus,         mPanPan1
    'GCA_029281585.2': 'gorilla',              # Gorilla gorilla,      mGorGor1
    'GCA_028885655.2': 'Sumatran_orangutan',   # Pongo abelii,         mPonAbe1
    'GCA_028885625.2': 'Bornean_orangutan',    # Pongo pygmaeus,       mPonPyg2
    'GCA_028878055.2': 'siamang',              # Symphalangus syndactylus, mSymSyn1
    'hs1': 'human_CHM13',
    'hg38': 'human_GRCh38',
}

# Guide tree used for the alignment, siamang as outgroup (README). Divergence
# order from human, for polarising an inversion onto a branch.
OUTGROUP_ORDER = ['chimpanzee', 'bonobo', 'gorilla',
                  'Sumatran_orangutan', 'Bornean_orangutan', 'siamang']

GCA_ASSEMBLIES = [a for a in SPECIES if a.startswith('GCA_')]

PRIMARY_CHROMS = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

# UCSC chain header: chain score tName tSize tStrand tStart tEnd
#                          qName qSize qStrand qStart qEnd id
CHAIN_COLS = ['kw', 'score', 'tName', 'tSize', 'tStrand', 'tStart', 'tEnd',
              'qName', 'qSize', 'qStrand', 'qStart', 'qEnd', 'id']


def chain_url(ref: str, query: str) -> str:
    """Chains are published as <target>_vs_<query>.chain.gz."""
    return f"{CHAIN_BASE}/{ref}_vs_{query}.chain.gz"


def chromalias_url(accession: str) -> str:
    acc = accession.split('.')[0]           # GCA_028858775
    digits = acc.split('_')[1]              # 028858775
    path = f"{digits[0:3]}/{digits[3:6]}/{digits[6:9]}"
    return (f"https://hgdownload.soe.ucsc.edu/hubs/GCA/{path}/"
            f"{accession}/{accession}.chromAlias.txt")
