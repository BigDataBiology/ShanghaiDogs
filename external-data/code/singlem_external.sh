#!/bin/bash

set -e
set -v

source ~/miniconda3/bin/activate root
conda activate kingfisher

cd /data/Projects/ShanghaiDogs/external-data/

# retrieve Biosample ID for each SRA run
kingfisher annotate --run-identifiers-list data/dog_microbiome_archive_otu_tables/SRR_Acc_List.txt \
--all-columns -o data/dog_microbiome_archive_otu_tables/SRR_Acc_List_Metadata.txt

conda deactivate
conda activate python_env

python code/link_run_to_biosample.py

conda deactivate
conda activate singleM

cd /data/Projects/ShanghaiDogs/external-data/data/dog_microbiome_archive_otu_tables/run_to_biosample
for ls in *.txt
  do
    prefix=$(echo "$ls" | cut -d '_' -f 1)
    head $ls
    echo $prefix
    singlem summarise --input-archive-otu-table $(cat $ls) \
    --output-archive-otu-table ../otu_tab_by_biosample/$prefix.json \
    --collapse-to-sample-name $prefix
  done

# Reannotate the archive OTU tables from Sandpiper
cd /data/Projects/ShanghaiDogs/external-data/data/dog_microbiome_archive_otu_tables/otu_tab_by_biosample
for json in *.json
  do
    prefix=$(echo "$json" | cut -d '.' -f 1)
    singlem renew --input-archive-otu-table ${json} \
    -p ${prefix}_profile.tsv \
    --taxonomic-profile-krona ${prefix}_krona.html
  done

