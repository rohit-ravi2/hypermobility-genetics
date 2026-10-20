#!/usr/bin/env bash
# Download the public GWAS summary statistics used here into data/sumstats/.
set -eu
cd "$(dirname "$0")/.."
mkdir -p data/sumstats && cd data/sumstats
dl(){ curl -sSL --retry 3 -o "$2" "$1"; echo "$2 $(stat -c%s "$2")"; }

# FinnGen R13 (hypermobility = ICD-10 M35.7)
F=https://storage.googleapis.com/finngen-public-data-r13/summary_stats
dl $F/finngen_R13_M13_HYPERMOBILITY.gz finngen_R13_hypermobility.gz
dl $F/finngen_R13_K11_HERNIA.gz        finngen_R13_hernia.gz
dl $F/finngen_R13_I9_VARICVE.gz        finngen_R13_varicose.gz
dl $F/finngen_R13_N14_FEMGENPROL.gz    finngen_R13_prolapse.gz

# PGC psychiatric GWAS (figshare)
dl https://ndownloader.figshare.com/files/40036684 ADHD2022_iPSYCH_deCODE_PGC.meta.gz
dl https://ndownloader.figshare.com/files/28169292 iPSYCH-PGC_ASD_Nov2017.gz
dl https://ndownloader.figshare.com/files/62089324 ANX_2026_daner_fullANX.gz
dl https://ndownloader.figshare.com/files/52419704 pgc-mdd2025_no23andMe_eur.tsv.gz

# Common executive function, Hatoum et al. 2023 (GWAS Catalog GCST90162547), LDSC sanity check
dl https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90162001-GCST90163000/GCST90162547/GCST90162547_buildGRCh37.tsv execfunction_hatoum2023.tsv
