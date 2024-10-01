#!/bin/bash
source ~/miniconda3/bin/activate root
conda activate prodigal

set -e
set -v

contigs="../data/ShanghaiDogsAssemblies"
output_folder="../intermediate-outputs/Prodigal"
temp_folder=$(mktemp -d)

for n in {00..52}
  do
    if [ "$n" -eq "09" ]; then
      continue
    fi
    echo D0${n}
    mkdir ${output_folder}/D0${n}/
    temp_fna=$temp_folder/D0${n}.fna
    zcat ${contigs}/D0${n}_PP1_PolcaCorr.fna.gz >$temp_fna
    prodigal -i ${temp_fna} \
        -o ${output_folder}/D0${n}/D0${n}_coords.gbk \
        -a ${output_folder}/D0${n}/D0${n}_proteins.faa \
        -d ${output_folder}/D0${n}/D0${n}_orfs.fna -p meta

    rm $temp_fna
  done
rm -rf $temp_folder
