#!/bin/bash

source ~/miniconda3/bin/activate root
conda activate masurca

out_PP1="/data/anna/animal_metagenome/long-mg-dog/02_polishing/01_polypolish"
out_path="/data/anna/animal_metagenome/long-mg-dog/02_polishing/02_POLCA"
input_reads="/data/anna/animal_metagenome/long-mg-dog/00_quality_control/01_trim_filter_data/ngs"

for n in {00..52}
  do
    mkdir ${out_path}/D0${n}
    cd ${out_path}/D0${n}

    polca.sh -a ${out_PP1}/D0${n}/D0${n}_PP_1.fasta \
    -r "${input_reads}/D0${n}_350_trim_filter.pair.1.fq.gz ${input_reads}/D0${n}_350_trim_filter.pair.2.fq.gz" \
    -t 32 -m 2G

    mv ${out_path}/D0${n}/D0${n}_PP_1.fasta.PolcaCorrected.fa ${out_path}/D0${n}/D0${n}_PP1_PolcaCorr.fasta
  done
