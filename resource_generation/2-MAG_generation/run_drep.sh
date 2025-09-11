#!/bin/bash

out_dastool="/data/anna/animal_metagenome/long-mg-dog/05_dereplication/00_dastool"
out_drep="/data/anna/animal_metagenome/long-mg-dog/05_dereplication/01_drep"

for n in {00..52}
  do
    suffix="_D0${n}"
    echo $suffix
    cd ${out_dastool}/D0${n}_DASTool_bins
    for file in *.fa.gz
      do
        filename=$(basename "$file" | cut -f 1 -d '.')
        # echo $filename
        extension=".fa.gz"
        # Rename the file with the sample_id
        mv "$file" "$filename$suffix$extension"
      done
  done

source ~/miniconda3/bin/activate root
conda activate dRep

# species-level dereplication
# conda update drep #it was not dereplicating before

cd ${out_drep}

dRep dereplicate -p 32 --genomeInfo quality_report_dastool_basename.csv \
--completeness 50 --contamination 10 ANI_95

dRep dereplicate -p 32 --genomeInfo quality_report_dastool_basename.csv \
--completeness 50 --contamination 10 --S_ani 0.995 ANI_995

dRep dereplicate -p 32 --genomeInfo quality_report_dastool_basename.csv \
--completeness 50 --contamination 10 --S_ani 0.9999 ANI_9999
