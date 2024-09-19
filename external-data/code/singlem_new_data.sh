#!/bin/bash

set -e
set -v

source ~/miniconda3/bin/activate root
conda activate singleM
export SINGLEM_METAPACKAGE_PATH='/data/yiqian/databases/singlem/S3.2.1.GTDB_r214.metapackage_20231006.smpkg.zb'

# Bioproject PRJNA871950
reads_path='/data/Projects/ShanghaiDogs/external-data/data/PRJNA871950/fastq'
cd /data/Projects/ShanghaiDogs/external-data/data/PRJNA871950

for fq in fastq/*.gz
  do
    a=$(basename $fq | sed 's/_.*//')
    mkdir singlem_profiles/"$a"
    singlem pipe -1 fastq/${a}_1.fastq.gz -2 fastq/${a}_2.fastq.gz \
    --threads 24 -p singlem_profiles/"$a"/"$a".profile.tsv \
    --otu-table singlem_profiles/"$a"/"$a"_otutable.tsv \
    --taxonomic-profile-krona singlem_profiles/"$a"/"$a"_krona.html
  done

# Bioproject
reads_path='/data/Projects/ShanghaiDogs/external-data/data/PRJNA917802/fastq'
cd /data/Projects/ShanghaiDogs/external-data/data/PRJNA917802

for fq in fastq/*.gz
  do
    a=$(basename $fq | sed 's/_.*//')
    mkdir singlem_profiles/"$a"
    singlem pipe -1 fastq/${a}_1.fastq.gz -2 fastq/${a}_2.fastq.gz \
    --threads 24 -p singlem_profiles/"$a"/"$a".profile.tsv \
    --otu-table singlem_profiles/"$a"/"$a"_otutable.tsv \
    --taxonomic-profile-krona singlem_profiles/"$a"/"$a"_krona.html
  done