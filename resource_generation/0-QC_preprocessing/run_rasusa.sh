#!/bin/bash

source ~/miniconda3/bin/activate root
conda activate rasusa

ont_filtered="/data/anna/animal_metagenome/long-mg-dog/00_quality_control/01_trim_filter_data/ont/porechop"
ont_subsets="/data/anna/animal_metagenome/long-mg-dog/00_quality_control/01_trim_filter_data/ont/porechop/subsets"

for n in {00..52}
  do
    rasusa --input ${ont_filtered}/D0${n}_filt_porechop.fastq.gz --bases 10G --output ${ont_subsets}/D0${n}_10G_filt_porechop.fastq.gz
  done

for n in {00..52}
  do
    rasusa --input ${ont_filtered}/D0${n}_filt_porechop.fastq.gz --bases 20G --output ${ont_subsets}/D0${n}_20G_filt_porechop.fastq.gz
  done
