#!/bin/bash

input_folder="/data/Projects/ShanghaiDogs/intermediate-outputs/Prodigal"
out_folder="/data/Projects/ShanghaiDogs/intermediate-outputs/eggNOG_annot_contigs"

source ~/miniconda3/bin/activate root
conda activate eggnog-mapper

export EGGNOG_DATA_DIR=/data/anna/miniconda3/envs/eggnog-mapper/db

for n in {00..52}
  do
    if [ "$n" -eq "09" ]; then
      continue
    fi
    echo D0${n}
    mkdir ${out_folder}/D0${n}/
    cd ${out_folder}/D0${n}/
    emapper.py -m diamond --itype proteins --cpu 32 \
    -i ${input_folder}/D0${n}/D0${n}_proteins.faa -o D0${n}
  done