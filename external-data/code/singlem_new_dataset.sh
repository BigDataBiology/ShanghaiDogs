#!/bin/bash

set -e
set -v

source ~/miniconda3/bin/activate root
conda activate singleM
export SINGLEM_METAPACKAGE_PATH=${HOME}/databases/singlem/S3.2.1.GTDB_r214.metapackage_20231006.smpkg.zb

# Bioproject PRJNA871950
reads_path='external-data/data/PRJNA871950'
cd $reads_path
output_path='intermediate-outputs/singlem_profiling/EXTERNAL-DATASETS/PRJNA871950'

for fq in fastq/*.gz
  do
    a=$(basename $fq | sed 's/_.*//')
    mkdir ${output_path}/"$a"
    singlem pipe -1 fastq/${a}_1.fastq.gz -2 fastq/${a}_2.fastq.gz \
    --threads 24 -p ${output_path}/"$a"/"$a".profile.tsv \
    --otu-table ${output_path}/"$a"/"$a"_otutable.tsv \
    --taxonomic-profile-krona ${output_path}/"$a"/"$a"_krona.html
  done

# Bioproject PRJNA917802
reads_path='external-data/data/PRJNA917802/'
cd ${reads_path}
output_path='intermediate-outputs/singlem_profiling/EXTERNAL-DATASETS/PRJNA917802'

for fq in fastq/*.gz
  do
    a=$(basename $fq | sed 's/_.*//')
    mkdir ${output_path}/"$a"
    singlem pipe -1 fastq/${a}_1.fastq.gz -2 fastq/${a}_2.fastq.gz \
    --threads 24 -p ${output_path}/"$a"/"$a".profile.tsv \
    --otu-table ${output_path}/"$a"/"$a"_otutable.tsv \
    --taxonomic-profile-krona ${output_path}/"$a"/"$a"_krona.html
  done

# Berlin dogs [unpublished data]
reads_path='external-data/data/berlin_dogs'
output_path='intermediate-outputs/singlem_profiling/EXTERNAL-DATASETS/berlin_dogs'
cd ${reads_path}

for fq in *.gz
  do
    a=$(basename $fq | sed 's/_.*//')
    file=$(basename $fq | sed 's/_[12].*//')
    mkdir ${output_path}/"$a"
    singlem pipe -1 ${file}_1.fq.gz -2 ${file}_2.fq.gz \
    --threads 24 -p ${output_path}/"$a"/"$a".profile.tsv \
    --otu-table ${output_path}/"$a"/"$a"_otutable.tsv \
    --taxonomic-profile-krona ${output_path}/"$a"/"$a"_krona.html
  done

# Bioproject NomNomNow (shallow metagenomics)
reads_path='external-data/data/USA_pets_Nomnomnow'
cd ${reads_path}
output_path='intermediate-outputs/singlem_profiling/EXTERNAL-DATASETS/USA_pets_Nomnomnow'

for fq in *.gz
  do
    a=$(basename $fq .fastq.gz)
    mkdir ${output_path}/"$a"
    singlem pipe -1 ${a}.fastq.gz \
    --threads 24 -p ${output_path}/"$a"/"$a".profile.tsv \
    --otu-table ${output_path}/"$a"/"$a"_otutable.tsv \
    --taxonomic-profile-krona ${output_path}/"$a"/"$a"_krona.html
  done

# Move post-treatment profiles into an independent folder (not included for SingleM)
mkdir ${output_path}/post_treatment_excluded/
mv SRR13972888 SRR13972889 SRR13972890 SRR13972891 SRR13972892 SRR13972893 SRR13972895 SRR13972896 SRR13972897 \
SRR13972898 SRR13972899 SRR13972900 SRR13972901 SRR13972902 SRR13972903 SRR13972904 SRR13972906 SRR13972907 \
SRR13972908 SRR13972909 SRR13972910 SRR13972911 SRR13972912 SRR13972913 SRR13972914 SRR13972915 SRR13972917 \
SRR13972918 ${output_path}/post_treatment_excluded/
