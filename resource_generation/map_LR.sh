#!/bin/bash

source ~/miniconda3/bin/activate root
conda activate minimap2

drep_path='/data/anna/animal_metagenome/long-mg-dog/05_dereplication/01_drep/ANI_95'

for n in {00..52}
  do
    minimap2 -t 32 -ax map-ont ${drep_path}/drep_95_MAGs_rf.fa \
    /data/anna/animal_metagenome/long-mg-dog/00_quality_control/01_trim_filter_data/ont/porechop/D0${n}_filt_porechop.fastq.gz |\
    samtools sort -o ${drep_path}/map_LR/D0${n}_LR_to_95_ANI.bam
    samtools flagstat ${drep_path}/map_LR/D0${n}_LR_to_95_ANI.bam
  done
