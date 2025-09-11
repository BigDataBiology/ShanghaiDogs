#!/bin/bash

source ~/miniconda3/bin/activate root
conda activate minimap2
input_path='/data/anna/animal_metagenome/long-mg-dog'
semibin_path='/data/anna/animal_metagenome/long-mg-dog/04_binning/00_SemiBin2/LR'

for n in {00..52}
  do
    mkdir ${semibin_path}/D0${n}
    minimap2 -ax map-ont ${input_path}/02_polishing/02_POLCA/D0${n}/D0${n}_PP1_PolcaCorr.fasta \
    ${input_path}/00_quality_control/01_trim_filter_data/ont/porechop/D0${n}_filt_porechop.fastq.gz > \
    ${semibin_path}/D0${n}/D0${n}_polca_long.sam -t 32
  done

conda deactivate
conda activate samtools

for n in {00..52}
  do
    samtools view -h -b -S ${semibin_path}/D0${n}/D0${n}_polca_long.sam -o ${semibin_path}/D0${n}/D0${n}_polca_long.bam
    samtools view -b -F 4 ${semibin_path}/D0${n}/D0${n}_polca_long.bam -o ${semibin_path}/D0${n}/D0${n}_polca_long_mapped.bam
    samtools sort ${semibin_path}/D0${n}/D0${n}_polca_long_mapped.bam -o ${semibin_path}/D0${n}/D0${n}_polca_long_mapped_sorted.bam
    samtools index ${semibin_path}/D0${n}/D0${n}_polca_long_mapped_sorted.bam
    samtools flagstat ${semibin_path}/D0${n}/D0${n}_polca_long.bam > ${semibin_path}/D0${n}/D0${n}_flagstat_long.txt
  done