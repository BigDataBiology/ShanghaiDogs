#!/bin/bash

source ~/miniconda3/bin/activate root
conda activate gunc

# SemiBin single SR
semibin_path='/data/anna/animal_metagenome/long-mg-dog/04_binning/00_SemiBin2'

for n in {00..52}
do
  echo "Starting D0'$n'..."
  cd "${semibin_path}/D0${n}/HQ_mq_MAGs" || continue
  for fa in *.fa.gz
  do
    out_dir=$(basename "$fa" .fa.gz)
    if [ ! -d "$out_dir" ]; then
      mkdir "$out_dir"
      gunc run -i "${fa}" --threads 32 \
      -r ~/miniconda3/envs/gunc/gunc_db/gunc_db_progenomes2.1.dmnd \
      --out_dir "$out_dir"
    else
      echo "Output directory '$out_dir' already exists. Skipping gunc."
    fi
  done
done

# SemiBin single LR
semibin_path='/data/anna/animal_metagenome/long-mg-dog/04_binning/00_SemiBin2/LR'

for n in {00..52}
do
  echo "Starting D0'$n'..."
  cd "${semibin_path}/D0${n}/HQ_mq_MAGs" || continue
  for fa in *.fa.gz
  do
    out_dir=$(basename "$fa" .fa.gz)
    if [ ! -d "$out_dir" ]; then
      mkdir "$out_dir"
      gunc run -i "${fa}" --threads 32 \
      -r ~/miniconda3/envs/gunc/gunc_db/gunc_db_progenomes2.1.dmnd \
      --out_dir "$out_dir"
    else
      echo "Output directory '$out_dir' already exists. Skipping gunc."
    fi
  done
done

# SemiBin single SR - subset 20G
semibin_path='/data/anna/animal_metagenome/long-mg-dog/04_binning/00_SemiBin2/subsets_SR'

for n in {00..52}
do
  echo "Starting D0'$n'..."
  cd "${semibin_path}/D0${n}_20G/HQ_mq_MAGs" || continue
  for fa in *.fa.gz
  do
    out_dir=$(basename "$fa" .fa.gz)
    if [ ! -d "$out_dir" ]; then
      mkdir "$out_dir"
      gunc run -i "${fa}" --threads 32 \
      -r ~/miniconda3/envs/gunc/gunc_db/gunc_db_progenomes2.1.dmnd \
      --out_dir "$out_dir"
    else
      echo "Output directory '$out_dir' already exists. Skipping gunc."
    fi
  done
done

# SemiBin single SR - subset 10G

for n in {00..52}
do
  echo "Starting D0'$n'..."
  cd "${semibin_path}/D0${n}_10G/HQ_mq_MAGs" || continue
  for fa in *.fa.gz
  do
    out_dir=$(basename "$fa" .fa.gz)
    if [ ! -d "$out_dir" ]; then
      mkdir "$out_dir"
      gunc run -i "${fa}" --threads 32 \
      -r ~/miniconda3/envs/gunc/gunc_db/gunc_db_progenomes2.1.dmnd \
      --out_dir "$out_dir"
    else
      echo "Output directory '$out_dir' already exists. Skipping gunc."
    fi
  done
done