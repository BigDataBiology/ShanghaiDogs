#!/bin/bash

ref="/data/anna/animal_metagenome/long-mg-dog/input/dog-ref-genome/GCA_000002285.4_Dog10K_Boxer_Tasha_genomic.fna"
resultpath="/data/anna/animal_metagenome/long-mg-dog/00_quality_control/01_trim_filter_data/ont/"
ont_path="/data/Public/Result-X101SC23043778-Z01-J005-B1-6_delivery_20230818/ONT"

samples_1="D000 D002 D003 D004 D005 D008"
ontdata_1="/data/anna/animal_metagenome/long-mg-dog/input/D000-D008/CP2019102500068/H101SC23043778/RSMD01614/X101SC23043778-Z01/X101SC23043778-Z01-F002/ONTdata_X101SC23043778-Z01-J001/"

samples_2="D001 D006"
ontdata_2="/data/anna/animal_metagenome/long-mg-dog/input/D000-D008/X101SC23043778-Z01/X101SC23043778-Z01-F001/ONTdata_X101SC23043778-Z01-J002"

samples_3="D007"
ontdata_3="/data/anna/animal_metagenome/long-mg-dog/input/D000-D008/wangxinyu0720067_D007_bis/Release-X101SC23043778-Z01-J004-20230715/Data-X101SC23043778-Z01-J004/D007/1506_6C_PAQ48914_ff9be250"

source ~/miniconda3/bin/activate root
conda activate chopper

for n in $samples_1;
  do
    gunzip -c ${ontdata_1}/${n}/*/*/${n}_pass.fastq.gz | chopper -q 10 -l 500 --threads 64 --contam ${ref} | gzip > ${resultpath}/${n}_ont_trim_filter.fastq.gz
  done

for n in $samples_2;
  do
    gunzip -c ${ontdata_2}/${n}/*/*/${n}_pass.fastq.gz | chopper -q 10 -l 500 --threads 64 --contam ${ref} | gzip > ${resultpath}/${n}_ont_trim_filter.fastq.gz
  done

for n in $samples_3;
  do
    gunzip -c ${ontdata_3}/${n}_pass.fastq.gz | chopper -q 10 -l 500 --threads 64 --contam ${ref} | gzip > ${resultpath}/${n}_ont_trim_filter.fastq.gz
  done

for n in {10..52}
  do
    gunzip -c ${ont_path}/D0${n}/D0${n}.fastq_pass.gz | chopper -q 10 -l 500 --threads 64 --contam ${ref} | gzip > ${resultpath}/D0${n}_ont_trim_filter.fastq.gz
  done
