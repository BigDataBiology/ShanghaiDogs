#!/bin/bash

source ~/miniconda3/bin/activate root
conda activate checkm
checkm data setRoot /data/anna/miniconda3/envs/checkm/checkm-db/

# Run Checkm1 to the MAG collection in the different polishing steps

input_path='/data/Projects/ShanghaiDogs/intermediate-outputs/polish_evaluation_mags'
output_path='/data/Projects/ShanghaiDogs/intermediate-outputs/polishing_evaluation/results/checkm1'

mkdir ${output_path}/flye
checkm lineage_wf -t 24 -x fna ${input_path}/flye_bins ${output_path}/flye

mkdir ${output_path}/medaka
checkm lineage_wf -t 24 -x fna ${input_path}/medaka_bins ${output_path}/medaka

mkdir ${output_path}/polypolish
checkm lineage_wf -t 24 -x fna ${input_path}/polypolish_bins ${output_path}/polypolish

mkdir ${output_path}/polca
checkm lineage_wf -t 24 -x gz /data/Projects/ShanghaiDogs/intermediate-outputs/polishing_evaluation/final_MAGs ${output_path}/polca

# Re-run Checkm1 to the bin collection (of all diff binning approaches)
input_path='/data/anna/animal_metagenome/long-mg-dog/04_binning'
output_path='/data/Projects/ShanghaiDogs/intermediate-outputs/polishing_evaluation/results/checkm1'

mkdir -p ${output_path}/SemiBin/ALL_LR_SemiBin
mkdir -p ${output_path}/SemiBin/ALL_SR_SemiBin
mkdir -p ${output_path}/SemiBin/subsets_SR_20G
mkdir -p ${output_path}/SemiBin/subsets_SR_10G
mkdir -p ${output_path}/SemiBin/ALL_multi_SemiBin

for n in {00..52}; do
  # Skip D009
  if [[ $n == "09" ]]; then
    continue
  fi
  mkdir -p ${output_path}/SemiBin/ALL_SR_SemiBin/D0${n}
  checkm lineage_wf -t 24 -x gz ${input_path}/00_SemiBin2/D0${n}/output_bins ${output_path}/SemiBin/ALL_SR_SemiBin/D0${n}
done

for n in {00..52}; do
  if [[ $n == "09" ]]; then
    continue
  fi
  mkdir ${output_path}/SemiBin/ALL_LR_SemiBin/D0${n}
  checkm lineage_wf -t 24 -x gz ${input_path}/00_SemiBin2/LR/D0${n}/output_bins ${output_path}/SemiBin/ALL_LR_SemiBin/D0${n}
done

for n in {00..52}; do
  if [[ $n == "09" ]]; then
    continue
  fi
  mkdir ${output_path}/SemiBin/subsets_SR_20G/D0${n}
  checkm lineage_wf -t 24 -x gz ${input_path}/00_SemiBin2/subsets_SR/D0${n}_20G/output_bins ${output_path}/SemiBin/subsets_SR_20G/D0${n}
done

for n in {00..52}; do
  if [[ $n == "09" ]]; then
    continue
  fi
  mkdir ${output_path}/SemiBin/subsets_SR_10G/D0${n}
  checkm lineage_wf -t 24 -x gz ${input_path}/00_SemiBin2/subsets_SR/D0${n}_10G/output_bins ${output_path}/SemiBin/subsets_SR_10G/D0${n}
done

for n in {00..52}; do
  if [[ $n == "09" ]]; then
    continue
  fi
  mkdir ${output_path}/SemiBin/ALL_multi_SemiBin/D0${n}
  checkm lineage_wf -t 24 -x gz ${input_path}/01_SemiBin2_multi/output_bins/D0${n}/output_bins ${output_path}/SemiBin/ALL_multi_SemiBin/D0${n}
done