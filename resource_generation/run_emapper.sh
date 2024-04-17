#!/bin/bash

input_folder="/data/Projects/ShanghaiDogs/data/ShanghaiDogsMAGs"
out_folder="/data/anna/animal_metagenome/long-mg-dog/09_eggNOG/eggNOG_annot"

source ~/miniconda3/bin/activate root
conda activate eggnog-mapper
export EGGNOG_DATA_DIR=/data/anna/miniconda3/envs/eggnog-mapper/db

cd ${input_folder}
for fafile in *.fna.gz
  do
    fbname=$(basename "$fafile" .fna.gz)
    echo $fbname
    emapper.py -m diamond --itype metagenome --genepred prodigal --cpu 32 \
    -i ${fafile} -o ${fbname} --output_dir ${out_folder}/
  done
