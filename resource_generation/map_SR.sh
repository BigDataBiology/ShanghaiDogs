#!/bin/bash

source ~/miniconda3/bin/activate root
conda activate strobealign

path="/data/anna/animal_metagenome/long-mg-dog"

for n in {00..52}
  do
    strobealign -t 8 ${path}/05_dereplication/01_drep/ANI_95/drep_95_MAGs_rf.fa \
    ${path}/00_quality_control/01_trim_filter_data/ngs/D0${n}_350_trim_filter.pair.1.fq.gz \
    ${path}/00_quality_control/01_trim_filter_data/ngs/D0${n}_350_trim_filter.pair.2.fq.gz \
    | samtools sort -o ${path}/05_dereplication/01_drep/ANI_95/align_SR/D0${n}_SR_to_95_ANI.bam
    samtools flagstat ${path}/05_dereplication/01_drep/ANI_95/align_SR/D0${n}_SR_to_95_ANI.bam
  done
