#!/bin/bash
source ~/miniconda3/bin/activate root
conda activate rgi

rgi load --card_json ~/miniconda3/envs/rgi/card.json --local

dna_contigs="/data/Projects/ShanghaiDogs/data/ShanghaiDogsAssemblies"
output_folder="/data/Projects/ShanghaiDogs/intermediate-outputs/06_ARG/01_RGI_CARD_contigs"

# DNA-level ARG prediction

for n in {00..52}
  do
    echo D0${n}
    mkdir ${output_folder}/D0${n}_dna/
    rgi main --input_sequence ${dna_contigs}/D0${n}_PP1_PolcaCorr.fna.gz \
    --output_file ${output_folder}/D0${n}_dna/out_D0${n} -t contig --clean
  done

for n in {00..52}
  do
    grep -h -e "Strict" -e "Perfect" ${output_folder}/D0${n}_dna/*.txt | \
    awk -v sample_id="D0${n}" '{sub(/\r$/, ""); print $0 "\t" sample_id}' > \
    ${output_folder}/D0${n}_dna/D0${n}_ARGs_formatted.txt
  done

cat ${output_folder}/D0*/*ARGs_formatted.txt > ${output_folder}/contigs-ARGs_ALL_samples.txt
