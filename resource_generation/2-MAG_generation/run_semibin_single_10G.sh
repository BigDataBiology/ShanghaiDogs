#!/bin/bash

source ~/miniconda3/bin/activate root
conda activate SemiBin
semibin_path='/data/anna/animal_metagenome/long-mg-dog/04_binning/00_SemiBin2/subsets_SR'

# Binning
for n in {00..52}
  do
    SemiBin2 single_easy_bin --environment dog_gut \
    -i /data/anna/animal_metagenome/long-mg-dog/02_polishing/02_POLCA/subsets/D0${n}_10G/D0${n}_PP1_PolcaCorr.fasta \
    -b ${semibin_path}/D0${n}_10G/D0${n}_polca_short_mapped_sorted.bam -o ${semibin_path}/D0${n}_10G \
    --sequencing-type long_read
  done

# Checkm2
conda deactivate
conda activate checkm2

for n in {00..52}
  do
    checkm2 predict -x gz --threads 24 --input ${semibin_path}/D0${n}_10G/output_bins \
    --output-directory ${semibin_path}/D0${n}_10G/checkm2/ --force
  done

# GTDB-tk
conda deactivate
conda activate gtdbtk

for n in {00..52}
  do
    gtdbtk classify_wf --cpus 24 -x gz --mash_db /data/anna/miniconda3/envs/gtdbtk/share/gtdbtk-2.3.0/db/mash/ \
    --genome_dir ${semibin_path}/D0${n}_10G/output_bins/ --out_dir ${semibin_path}/D0${n}_10G/gtdbtk/
  done
