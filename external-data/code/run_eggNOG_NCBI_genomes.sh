#!/bin/bash

set -e
set -v

input_folder="/data/Projects/ShanghaiDogs/external-data/data/NCBI_genomes_ref"
out_folder="/data/Projects/ShanghaiDogs/external-data/data/NCBI_genomes_ref/eggNOG_annot"

source ~/miniconda3/bin/activate root
conda activate eggnog-mapper
export EGGNOG_DATA_DIR=/data/anna/miniconda3/envs/eggnog-mapper/db

cd ${input_folder}/refseq/bacteria/
for fafile in */*.fna.gz
  do
    fbname=$(basename "$fafile" .fna.gz)
    echo $fbname
    emapper.py -m diamond --itype metagenome --genepred prodigal --cpu 32 \
    -i ${fafile} -o ${fbname} --output_dir ${out_folder}/
  done

cd ${input_folder}/genbank/bacteria/
for fafile in */*.fna.gz
  do
    fbname=$(basename "$fafile" .fna.gz)
    echo $fbname
    emapper.py -m diamond --itype metagenome --genepred prodigal --cpu 32 \
    -i ${fafile} -o ${fbname} --output_dir ${out_folder}/
  done
