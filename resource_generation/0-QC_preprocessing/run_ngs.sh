#!/bin/bash

source ~/miniconda3/bin/activate root
conda activate ngless

ngs_data_0="/data/Public/X101SC23043778-Z01/Clean_Illumina_AUG_2023/X101SC23043778-Z01-F008"
ngs_data_1="/data/Public/Result-X101SC23043778-Z01-J005-B1-6_delivery_20230818/NGS/01.CleanData"
result_path="/data/anna/animal_metagenome/long-mg-dog/00_quality_control/01_trim_filter_data/ngs"

for n in {00..08}
  do
    ngless -j 64 trim_filter_ngs.ngl  ${ngs_data_0}/D0${n}/D0${n}_350.fq1.gz \
    ${ngs_data_0}/D0${n}/D0${n}_350.fq2.gz ${result_path}/D0${n}_350
  done

for n in {10..52}
  do
    ngless -j 64 trim_filter_ngs.ngl ${ngs_data_1}/D0${n}/D0${n}_350.fq1.gz \
    ${ngs_data_1}/D0${n}/D0${n}_350.fq2.gz ${result_path}/D0${n}_350
  done
