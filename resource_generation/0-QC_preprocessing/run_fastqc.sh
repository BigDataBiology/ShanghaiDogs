#!/bin/bash

source ~/miniconda3/bin/activate root
conda activate fastqc

ngs_data_1='/data/Public/X101SC23043778-Z01/Clean_Illumina_AUG_2023/X101SC23043778-Z01-F008'
ngs_data_2='/data/Public/Result-X101SC23043778-Z01-J005-B1-6_delivery_20230818/NGS/01.CleanData'
ngs_results='/data/anna/animal_metagenome/long-mg-dog/00_quality_control/00_raw_clean_data/ngs/00_fastQC'

ont_data_1='/data/Public/ANNA/wangxinyu0616030/Release-X101SC23043778-Z01-J001-B1-6_20230615/ONTdata_X101SC23043778-Z01-J001'
ont_data_2='/data/Public/Result-X101SC23043778-Z01-J005-B1-6_delivery_20230818/ONT'
ont_results='/data/anna/animal_metagenome/long-mg-dog/00_quality_control/00_raw_clean_data/ont/00_fastQC'

### NGS
for n in {00..08}
  do
    fastqc -o ${ngs_results} ${ngs_data_1}/D0${n}/D0${n}_350.fq1.gz ${ngs_data}/D0${n}/D0${n}_350.fq2.gz -t 64
  done

for n in {10..52}
  do
    fastqc -o ${ngs_results} ${ngs_data_2}/D0${n}/D0${n}_350.fq1.gz ${ngs_data}/D0${n}/D0${n}_350.fq2.gz -t 64
  done

### ONT

for n in {00..08}
  do
    fastqc -o ${ont_results} ${ont_data_1}/D0${n}/*/*/D0${n}.fastq_pass.gz -t 64
  done

# D001
fastqc -o ${ont_results} /data/Public/ANNA/wangxinyu0616030/Release-X101SC23043778-Z01-J002-B1-6_20230619/\
ONTdata_X101SC23043778-Z01-J002/D001/*/*/D001.fastq_pass.gz -t 64

# D006
fastqc -o ${ont_results} /data/Public/ANNA/wangxinyu0616030/Release-X101SC23043778-Z01-J002-B1-6_20230619/\
ONTdata_X101SC23043778-Z01-J002/D006/*/*/D006.fastq_pass.gz -t 64

# D007 (sample repetition)
fastqc -o ${ont_results} /data/Public/X101SC23043778-Z01/D007_new/D007.fastq_pass.gz -t 64

# New samples
for n in {10..52}
  do 
    fastqc -o ${ont_results} ${ont_data_2}/D0${n}/D0${n}.fastq_pass.gz -t 64
  done
