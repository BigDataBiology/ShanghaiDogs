#!/bin/bash
source ~/miniconda3/bin/activate root
conda activate rgi

rgi load --card_json ~/miniconda3/envs/rgi/card.json --local
out_dastool="/data/anna/animal_metagenome/long-mg-dog/05_dereplication/00_dastool"

for n in {00..52}
  do
    gzip -d ${out_dastool}/D0${n}_DASTool_bins/*.fa.gz
  done

for n in {00..52}
  do
    mkdir ~/animal_metagenome/long-mg-dog/06_ARG/00_RGI_CARD/D0${n}/
    cd ${out_dastool}/D0${n}_DASTool_bins
    for bin in *.fa
      do
        echo D0${n}
        echo ${bin}
        rgi main --input_sequence ${bin} --output_file ~/animal_metagenome/long-mg-dog/06_ARG/00_RGI_CARD/D0${n}/out_${bin} --clean --include_loose
      done
  done

for n in {00..52}
  do
    grep -e "Strict" -e "Perfect" /data/anna/animal_metagenome/long-mg-dog/06_ARG/00_RGI_CARD/D0${n}/*.txt > \
    /data/anna/animal_metagenome/long-mg-dog/06_ARG/00_RGI_CARD/D0${n}/D0${n}_ARG_strict_matches.txt
  done

for n in {00..52}
  do
    gzip ${out_dastool}/D0${n}_DASTool_bins/*.fa
  done