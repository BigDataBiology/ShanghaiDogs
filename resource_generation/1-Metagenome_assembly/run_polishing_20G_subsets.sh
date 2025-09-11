#!/bin/bash

source ~/miniconda3/bin/activate root

ont_reads="/data/anna/animal_metagenome/long-mg-dog/00_quality_control/01_trim_filter_data/ont/porechop/subsets"
flye_assembly="/data/anna/animal_metagenome/long-mg-dog/01_assembly/subset_assemblies"
out_medaka="/data/anna/animal_metagenome/long-mg-dog/02_polishing/00_medaka/subsets"
ngs_reads="/data/anna/animal_metagenome/long-mg-dog/00_quality_control/01_trim_filter_data/ngs"
out_polypolish="/data/anna/animal_metagenome/long-mg-dog/02_polishing/01_polypolish/subsets"
out_polca="/data/anna/animal_metagenome/long-mg-dog/02_polishing/02_POLCA/subsets"

# medaka
conda activate medaka-cpu

for n in {00..06}
  do
    medaka_consensus -i ${ont_reads}/D0${n}_20G_filt_porechop.fastq.gz \
    -d ${flye_assembly}/D0${n}_20G/assembly.fasta -o ${out_medaka}/D0${n}_20G/ -t 8 -m r1041_e82_400bps_hac_g632
  done

n=08
medaka_consensus -i ${ont_reads}/D0${n}_20G_filt_porechop.fastq.gz \
-d ${flye_assembly}/D0${n}_20G/assembly.fasta -o ${out_medaka}/D0${n}_20G/ -t 8 -m r1041_e82_400bps_hac_g632

# D007 and from D010 to D052 were basecalled with a newer Guppy version
n=07
medaka_consensus -i ${ont_reads}/D0${n}_20G_filt_porechop.fastq.gz \
-d ${flye_assembly}/D0${n}_20G/assembly.fasta -o ${out_medaka}/D0${n}_20G/ -t 8 -m r1041_e82_400bps_hac_v4.2.0

for n in {10..52}
  do
    medaka_consensus -i ${ont_reads}/D0${n}_20G_filt_porechop.fastq.gz \
    -d ${flye_assembly}/D0${n}_20G/assembly.fasta -o ${out_medaka}/D0${n}_20G/ -t 8 -m r1041_e82_400bps_hac_v4.2.0
   done

# Polypolish

conda deactivate
conda activate bwa

for n in {00..52}
  do
    bwa index ${out_medaka}/D0${n}_20G/consensus.fasta

    bwa mem -t 32 -a ${out_medaka}/D0${n}_20G/consensus.fasta ${ngs_reads}/D0${n}_350_trim_filter.pair.1.fq.gz \
    > ${out_polypolish}/D0${n}_20G/D0${n}_alignments_1.sam
    bwa mem -t 32 -a ${out_medaka}/D0${n}_20G/consensus.fasta ${ngs_reads}/D0${n}_350_trim_filter.pair.2.fq.gz > \
    ${out_polypolish}/D0${n}_20G/D0${n}_alignments_2.sam

    mkdir ${out_polypolish}/D0${n}_20G/
  done

conda deactivate
conda activate polypolish

for n in {00..52}
  do
    polypolish_insert_filter.py --in1 ${out_polypolish}/D0${n}_20G/D0${n}_alignments_1.sam \
    --in2 ${out_polypolish}/D0${n}_20G/D0${n}_alignments_2.sam \
    --out1 ${out_polypolish}/D0${n}_20G/D0${n}_filtered_1.sam \
    --out2 ${out_polypolish}/D0${n}_20G/D0${n}_filtered_2.sam
    polypolish ${out_medaka}/D0${n}_20G/consensus.fasta \
    ${out_polypolish}/D0${n}_20G/D0${n}_filtered_1.sam \
    ${out_polypolish}/D0${n}_20G/D0${n}_filtered_2.sam \
    > ${out_polypolish}/D0${n}_20G/D0${n}_PP_1.fasta
  done

# POLCA

conda deactivate
conda activate masurca

for n in {00..52}
  do
    mkdir ${out_polca}/D0${n}_20G
    cd ${out_polca}/D0${n}_20G

    polca.sh -a ${out_polypolish}/D0${n}_20G/D0${n}_PP_1.fasta \
    -r "${ngs_reads}/D0${n}_350_trim_filter.pair.1.fq.gz ${ngs_reads}/D0${n}_350_trim_filter.pair.2.fq.gz" \
    -t 32 -m 2G

    mv ${out_polca}/D0${n}_20G/D0${n}_PP_1.fasta.PolcaCorrected.fa ${out_path}/D0${n}_20G/D0${n}_PP1_PolcaCorr.fasta
  done

