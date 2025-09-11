#!/bin/bash

source ~/miniconda3/bin/activate root
conda activate python_env

ont_data="/data/anna/animal_metagenome/long-mg-dog/00_quality_control/01_trim_filter_data/ont"
result_path="/data/anna/animal_metagenome/long-mg-dog/00_quality_control/01_trim_filter_data/ont/porechop"

for n in {00..52}
  do
    porechop_abi -abi -i ${ont_data}/D0${n}_ont_trim_filter.fastq.gz --format fastq.gz --verbosity 1 \
    --discard_middle -t 32 -o ${result_path}/D0${n}_filt_porechop.fastq.gz > ${result_path}/out_D0${n}.txt
  done

