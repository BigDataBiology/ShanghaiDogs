#!/bin/bash

source ~/miniconda3/bin/activate root
conda activate medaka-cpu

ont_reads="/data/anna/animal_metagenome/long-mg-dog/00_quality_control/01_trim_filter_data/ont/porechop"
flye_assembly="/data/anna/animal_metagenome/long-mg-dog/01_assembly"
out_medaka="/data/anna/animal_metagenome/long-mg-dog/02_polishing/00_medaka"

for n in {00..06}
  do
    medaka_consensus -i ${ont_reads}/D0${n}_filt_porechop.fastq.gz \
    -d ${flye_assembly}/D0${n}/assembly.fasta -o ${out_medaka}/D0${n}_bo/ -t 6 -m r1041_e82_400bps_hac_g632
  done

n=08
medaka_consensus -i ${ont_reads}/D0${n}_filt_porechop.fastq.gz \
-d ${flye_assembly}/D0${n}/assembly.fasta -o ${out_medaka}/D0${n}_bo/ -t 6 -m r1041_e82_400bps_hac_g632


# D007 and from D010 to D052 were basecalled with a newer Guppy version
n=07
medaka_consensus -i ${ont_reads}/D0${n}_filt_porechop.fastq.gz \
-d ${flye_assembly}/D0${n}/assembly.fasta -o ${out_medaka}/D0${n}_bo/ -t 6 -m r1041_e82_400bps_hac_v4.2.0

for n in {10..52}
  do
    medaka_consensus -i ${ont_reads}/D0${n}_filt_porechop.fastq.gz \
    -d ${flye_assembly}/D0${n}/assembly.fasta -o ${out_medaka}/D0${n}_bo/ -t 6 -m r1041_e82_400bps_hac_v4.2.0
  done
