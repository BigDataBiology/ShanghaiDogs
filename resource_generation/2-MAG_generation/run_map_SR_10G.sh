#!/bin/bash

source ~/miniconda3/bin/activate root
conda activate bwa

input_path='/data/anna/animal_metagenome/long-mg-dog'
semibin_path='/data/anna/animal_metagenome/long-mg-dog/04_binning/00_SemiBin2/subsets_SR'

for n in {00..52}
  do
    bwa index ${input_path}/02_polishing/02_POLCA/subsets/D0${n}_10G/D0${n}_PP1_PolcaCorr.fasta
    bwa mem ${input_path}/02_polishing/02_POLCA/subsets/D0${n}_10G/D0${n}_PP1_PolcaCorr.fasta \
    ${input_path}/00_quality_control/01_trim_filter_data/ngs/D0${n}_350_trim_filter.pair.1.fq.gz \
    ${input_path}/00_quality_control/01_trim_filter_data/ngs/D0${n}_350_trim_filter.pair.2.fq.gz > \
    ${semibin_path}/D0${n}_10G/D0${n}_polca_short.sam -t 24
  done

conda deactivate
conda activate samtools

for n in {00..52}
  do
    samtools view -h -b -S ${semibin_path}/D0${n}_10G/D0${n}_polca_short.sam -o ${semibin_path}/D0${n}_10G/D0${n}_polca_short.bam
    samtools view -b -F 4 ${semibin_path}/D0${n}_10G/D0${n}_polca_short.bam -o ${semibin_path}/D0${n}_10G/D0${n}_polca_short_mapped.bam
    samtools sort ${semibin_path}/D0${n}_10G/D0${n}_polca_short_mapped.bam -o ${semibin_path}/D0${n}_10G/D0${n}_polca_short_mapped_sorted.bam
    samtools index ${semibin_path}/D0${n}_10G/D0${n}_polca_short_mapped_sorted.bam
    samtools flagstat ${semibin_path}/D0${n}_10G/D0${n}_polca_short.bam > ${semibin_path}/D0${n}_10G/D0${n}_flagstat_short.txt
  done