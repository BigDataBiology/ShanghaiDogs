#!/bin/bash

source ~/miniconda3/bin/activate root
conda activate Flye-2.9.2

ont_filtered="/data/anna/animal_metagenome/long-mg-dog/00_quality_control/01_trim_filter_data/ont/porechop"
out_flye="/data/anna/animal_metagenome/long-mg-dog/01_assembly"

# Error rate varies from 0.03-0.05. I will adjust accordingly. Assemblies should not be affected by 1-2% change in thresholds.
# Stated by developer in this Github issue: https://github.com/fenderglass/Flye/issues/508

declare -a list_03=("07" "10" "11" "12" "13" "14" "15" "16" "17" "19" "21" "22" "25" "27" "29" "30" "32" "34" "35" "38" "39" "40" "42" "44" "45" "46" "47" "48" "50" "51")
declare -a list_04=("00" "01" "03" "06" "08" "18" "20" "23" "26" "31" "36" "37" "41" "43" "49" "52" )
declare -a list_05=("02" "04" "05" "24" "28" "33" "46")

for n in ${list_03[@]}
  do
    echo ${n}
    flye -o ${out_flye}/D0${n}/ --threads 32 --nano-hq ${ont_filtered}/D0${n}_filt_porechop.fastq.gz --meta --read-error 0.03
  done

for n in ${list_04[@]}
  do
    echo ${n}
    flye -o ${out_flye}/D0${n}/ --threads 32 --nano-hq ${ont_filtered}/D0${n}_filt_porechop.fastq.gz --meta --read-error 0.04
  done

for n in ${list_05[@]}
  do
    echo ${n}
    flye -o ${out_flye}/D0${n}/ --threads 32 --nano-hq ${ont_filtered}/D0${n}_filt_porechop.fastq.gz --meta --read-error 0.05
  done