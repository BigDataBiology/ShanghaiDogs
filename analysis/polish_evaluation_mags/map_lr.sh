#!/bin/bash

for i in {0..8}
  do
    coverm genome --single /data/anna/animal_metagenome/long-mg-dog/00_quality_control/01_trim_filter_data/ont/porechop/D00${i}_filt_porechop.fastq.gz --genome-fasta-directory /data/Projects/ShanghaiDogs/intermediate-outputs/MAGs_by_sample/D00${i} --genome-fasta-extension gz --output-file /data/Projects/ShanghaiDogs/intermediate-outputs/polish_evaluation_mags/map_lr/D00${i}.tsv --output-format dense --mapper minimap2-ont --methods mean -t 40
  done

for i in {10..52}
  do 
    coverm genome --single /data/anna/animal_metagenome/long-mg-dog/00_quality_control/01_trim_filter_data/ont/porechop/D0${i}_filt_porechop.fastq.gz --genome-fasta-directory /data/Projects/ShanghaiDogs/intermediate-outputs/MAGs_by_sample/D0${i} --genome-fasta-extension gz --output-file /data/Projects/ShanghaiDogs/intermediate-outputs/polish_evaluation_mags/map_lr/D0${i}.tsv --output-format dense --mapper minimap2-ont --methods mean -t 40
  done
