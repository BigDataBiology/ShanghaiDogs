#!/bin/bash
source ~/miniconda3/bin/activate root
conda activate prodigal

set -e
set -v

contigs="/data/Projects/ShanghaiDogs/data/ShanghaiDogsAssemblies"
output_folder="/data/Projects/ShanghaiDogs/intermediate-outputs/Prodigal"

for n in {00..52}
  do
    if [ "$n" -eq "09" ]; then
      continue
    fi
    echo D0${n}
    mkdir ${output_folder}/D0${n}/
    gunzip ${contigs}/D0${n}_PP1_PolcaCorr.fna.gz
    prodigal -i ${contigs}/D0${n}_PP1_PolcaCorr.fna \
    -o ${output_folder}/D0${n}/D0${n}_coords.gbk \
    -a ${output_folder}/D0${n}/D0${n}_proteins.faa -p meta
    gzip ${contigs}/D0${n}_PP1_PolcaCorr.fna
  done