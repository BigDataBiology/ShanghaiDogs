#!/bin/bash

set -e
set -v

for n in {00..52}
  do
    ln -s /data/anna/animal_metagenome/long-mg-dog/05_dereplication/00_dastool/D0${n}_DASTool_bins/checkm2/quality_report.tsv \
    /data/Projects/ShanghaiDogs/intermediate-outputs/05_dereplication/00_dastool/D0${n}/checkm2/quality_report.tsv
    ln -s /data/anna/animal_metagenome/long-mg-dog/05_dereplication/00_dastool/D0${n}_DASTool_bins/gtdbtk/classify/gtdbtk.bac120.summary.tsv \
    /data/Projects/ShanghaiDogs/intermediate-outputs/05_dereplication/00_dastool/D0${n}/gtdbtk/gtdbtk.bac120.summary.tsv
    ln -s /data/anna/animal_metagenome/long-mg-dog/05_dereplication/00_dastool/D0${n}_DASTool_bins/subsets/checkm2/quality_report.tsv \
    /data/Projects/ShanghaiDogs/intermediate-outputs/05_dereplication/00_dastool/D0${n}/subsets/checkm2/quality_report.tsv
    ln -s /data/anna/animal_metagenome/long-mg-dog/05_dereplication/00_dastool/D0${n}_DASTool_bins/tigs_num_D0${n}.txt \
    /data/Projects/ShanghaiDogs/intermediate-outputs/05_dereplication/00_dastool/D0${n}/tigs_num_D0${n}.txt
  done

for n in {00..52}
  do
    mkdir /data/Projects/ShanghaiDogs/intermediate-outputs/05_dereplication/00_dastool/D0${n}/barrnap_out
    ln -s /data/anna/animal_metagenome/long-mg-dog/07_ribosomal_genes/barrnap_out/D0${n}/D0${n}_ALL.txt \
    /data/Projects/ShanghaiDogs/intermediate-outputs/05_dereplication/00_dastool/D0${n}/barrnap_out/D0${n}_ALL.txt
    mkdir /data/Projects/ShanghaiDogs/intermediate-outputs/05_dereplication/00_dastool/D0${n}/RGI_CARD
    ln -s /data/anna/animal_metagenome/long-mg-dog/06_ARG/00_RGI_CARD/D0${n}/D0${n}_ARG_strict_matches.txt.gz \
    /data/Projects/ShanghaiDogs/intermediate-outputs/05_dereplication/00_dastool/D0${n}/RGI_CARD/D0${n}_ARG_strict_matches.txt.gz
  done
