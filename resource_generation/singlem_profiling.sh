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
    # Calculate microbial fraction, size file created with singlem_metagenome_size.py python script
    singlem microbial_fraction --input-metagenome-sizes /data/anna/animal_metagenome/long-mg-dog/singlem_profiling/00_profiles/D0${n}/D0${n}_350_size.tsv \
    -p /data/anna/animal_metagenome/long-mg-dog/singlem_profiling/00_profiles/D0${n}/D0${n}.profile.tsv \
    --output-tsv /data/anna/animal_metagenome/long-mg-dog/singlem_profiling/00_profiles/D0${n}/D0${n}_microbial_fraction.tsv
  done

for n in {10..52}
  do
    mkdir /data/anna/animal_metagenome/long-mg-dog/singlem_profiling/00_profiles/D0${n}/
    singlem pipe -1 ${ngs_path_2}/D0${n}/D0${n}_350.fq1.gz \
    -2 ${ngs_path_2}/D0${n}/D0${n}_350.fq2.gz --threads 24 \
    -p /data/anna/animal_metagenome/long-mg-dog/singlem_profiling/00_profiles/D0${n}/D0${n}.profile.tsv \
    --otu-table /data/anna/animal_metagenome/long-mg-dog/singlem_profiling/00_profiles/D0${n}/D0${n}_otutable.tsv \
    --taxonomic-profile-krona /data/anna/animal_metagenome/long-mg-dog/singlem_profiling/00_profiles/D0${n}/D0${n}_krona.html
    # Calculate microbial fraction, size file created with singlem_metagenome_size.py python script
    singlem microbial_fraction --input-metagenome-sizes /data/anna/animal_metagenome/long-mg-dog/singlem_profiling/00_profiles/D0${n}/D0${n}_350_size.tsv \
    -p /data/anna/animal_metagenome/long-mg-dog/singlem_profiling/00_profiles/D0${n}/D0${n}.profile.tsv \
    --output-tsv /data/anna/animal_metagenome/long-mg-dog/singlem_profiling/00_profiles/D0${n}/D0${n}_microbial_fraction.tsv
  done

# Summarize taxonomic profiles at genus and species level

cd /data/Projects/ShanghaiDogs/

singlem summarise --input-taxonomic-profiles external-data/data/dog_microbiome_archive_otu_tables/otu_tab_by_biosample/done/*profile.tsv \
intermediate-outputs/singlem_profiling/D*/*profile.tsv \
--output-species-by-site-relative-abundance intermediate-outputs/singlem_profiling/all-dog-tax-profile-genus.tsv \
--output-species-by-site-level genus

singlem summarise --input-taxonomic-profiles external-data/data/dog_microbiome_archive_otu_tables/otu_tab_by_biosample/done/*profile.tsv \
intermediate-outputs/singlem_profiling/D*/*profile.tsv \
--output-species-by-site-relative-abundance intermediate-outputs/singlem_profiling/all-dog-tax-profile-species.tsv \
--output-species-by-site-level species

# Summarize otu tables

singlem summarise --input-archive-otu-tables external-data/data/dog_microbiome_archive_otu_tables/otu_tab_by_biosample/done/*.json \
--input-otu-tables intermediate-outputs/singlem_profiling/D*/*_otutable.tsv \
--output-otu-table intermediate-outputs/singlem_profiling/all-dog-combined-otu_table.csv


