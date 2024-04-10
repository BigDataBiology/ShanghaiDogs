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
  done


