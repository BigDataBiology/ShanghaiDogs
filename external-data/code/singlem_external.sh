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
    echo $prefix
    singlem summarise --input-archive-otu-table-list $ls \
    --output-archive-otu-table ../otu_tab_by_biosample/$prefix.json \
    --collapse-to-sample-name $prefix
  done
