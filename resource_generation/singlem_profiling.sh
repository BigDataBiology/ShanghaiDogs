#!/bin/bash

source ~/miniconda3/bin/activate root
conda activate singleM

export SINGLEM_METAPACKAGE_PATH='/data/yiqian/databases/singlem/S3.2.1.GTDB_r214.metapackage_20231006.smpkg.zb'
ngs_path_1='/data/Public/X101SC23043778-Z01/Clean_Illumina_AUG_2023/X101SC23043778-Z01-F008'
ngs_path_2='/data/Public/Result-X101SC23043778-Z01-J005-B1-6_delivery_20230818/NGS/01.CleanData'

for n in {00..08}
  do
    mkdir /data/anna/animal_metagenome/long-mg-dog/singlem_profiling/00_profiles/D0${n}/
    singlem pipe -1 ${ngs_path_1}/D0${n}/D0${n}_350.fq1.gz \
    -2 ${ngs_path_1}/D0${n}/D0${n}_350.fq2.gz --threads 24 \
    -p /data/anna/animal_metagenome/long-mg-dog/singlem_profiling/00_profiles/D0${n}/D0${n}.profile.tsv \
    --otu-table /data/anna/animal_metagenome/long-mg-dog/singlem_profiling/00_profiles/D0${n}/D0${n}_otutable.tsv \
    --taxonomic-profile-krona /data/anna/animal_metagenome/long-mg-dog/singlem_profiling/00_profiles/D0${n}/D0${n}_krona.html
  done

for n in {10..52}
  do
    mkdir /data/anna/animal_metagenome/long-mg-dog/singlem_profiling/00_profiles/D0${n}/
    singlem pipe -1 ${ngs_path_2}/D0${n}/D0${n}_350.fq1.gz \
    -2 ${ngs_path_2}/D0${n}/D0${n}_350.fq2.gz --threads 24 \
    -p /data/anna/animal_metagenome/long-mg-dog/singlem_profiling/00_profiles/D0${n}/D0${n}.profile.tsv \
    --otu-table /data/anna/animal_metagenome/long-mg-dog/singlem_profiling/00_profiles/D0${n}/D0${n}_otutable.tsv \
    --taxonomic-profile-krona /data/anna/animal_metagenome/long-mg-dog/singlem_profiling/00_profiles/D0${n}/D0${n}_krona.html
  done

# Summarize taxonomic profiles at genus and species level

cd /data/anna/animal_metagenome/long-mg-dog/singlem_profiling/00_profiles/

singlem summarise --input-taxonomic-profiles */*profile.tsv \
--output-species-by-site-relative-abundance tax-profile-genus.tsv \
--output-species-by-site-level genus

singlem summarise --input-taxonomic-profiles */*profile.tsv \
--output-species-by-site-relative-abundance tax-profile-species.tsv \
--output-species-by-site-level species

# Summarize otu tables

singlem summarise --input-otu-tables */*_otutable.tsv \
--output-otu-table combined.otu_table.csv