out_medaka="/data/anna/animal_metagenome/long-mg-dog/02_polishing/00_medaka"
input_reads="/data/anna/animal_metagenome/long-mg-dog/00_quality_control/01_trim_filter_data/ngs"
out_path="/data/anna/animal_metagenome/long-mg-dog/02_polishing/01_polypolish"

source ~/miniconda3/bin/activate root

conda activate bwa

for n in {00..52}
  do
    bwa index ${out_medaka}/D0${n}_bo/consensus.fasta
    bwa mem -t 32 -a ${out_medaka}/D0${n}_bo/consensus.fasta ${input_reads}/D0${n}_350_trim_filter.pair.1.fq.gz > \
    ${out_path}/D0${n}/D0${n}_alignments_1.sam
    bwa mem -t 32 -a ${out_medaka}/D0${n}_bo/consensus.fasta ${input_reads}/D0${n}_350_trim_filter.pair.2.fq.gz > \
    ${out_path}/D0${n}/D0${n}_alignments_2.sam
  done

conda deactivate
conda activate polypolish

for n in {00..52}
  do
    mkdir ${out_path}/D0${n}

    polypolish_insert_filter.py --in1 ${out_path}/D0${n}/D0${n}_alignments_1.sam \
    --in2 ${out_path}/D0${n}/D0${n}_alignments_2.sam \
    --out1 ${out_path}/D0${n}/D0${n}_filtered_1.sam \
    --out2 ${out_path}/D0${n}/D0${n}_filtered_2.sam

    polypolish ${out_medaka}/D0${n}_bo/consensus.fasta \
    ${out_path}/D0${n}/D0${n}_filtered_1.sam ${out_path}/D0${n}/D0${n}_filtered_2.sam > \
    ${out_path}/D0${n}/D0${n}_PP_1.fasta
  done
