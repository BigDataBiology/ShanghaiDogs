#!/bin/bash

source ~/miniconda3/bin/activate root
conda activate das_tool

binning_path='/data/anna/animal_metagenome/long-mg-dog/04_binning'

# Use the high and medium-quality MAGs that passed the GUNC filter as the QC_datasets
# Run on the MAGs that come from the same metagenome assembly and different binning (SR, LR, and multi)

for n in {00..52}
  do
    cd ${binning_path}/00_SemiBin2/D0${n}/HQ_mq_MAGs
    gzip -d *.fa.gz
    /data/anna/miniconda3/envs/das_tool/share/das_tool-1.1.6-0/src/Fasta_to_Contig2Bin.sh -e fa > contig_bins_dastool.tsv
    gzip *.fa
    cd ${binning_path}/00_SemiBin2/LR/D0${n}/HQ_mq_MAGs
    gzip -d *.fa.gz
    /data/anna/miniconda3/envs/das_tool/share/das_tool-1.1.6-0/src/Fasta_to_Contig2Bin.sh -e fa > contig_bins_dastool.tsv
    gzip *.fa
    cd ${binning_path}/01_SemiBin2_multi/output_bins/D0${n}/HQ_mq_MAGs
    gzip -d *.fa.gz
    /data/anna/miniconda3/envs/das_tool/share/das_tool-1.1.6-0/src/Fasta_to_Contig2Bin.sh -e fa > contig_bins_dastool.tsv
    gzip *.fa
  done

for n in {00..52}
  do
    DAS_Tool -i ${binning_path}/00_SemiBin2/D0${n}/HQ_mq_MAGs/contig_bins_dastool.tsv,${binning_path}/00_SemiBin2/LR/D0${n}/HQ_mq_MAGs/contig_bins_dastool.tsv,${binning_path}/01_SemiBin2_multi/output_bins/D0${n}/HQ_mq_MAGs/contig_bins_dastool.tsv \
    -c /data/anna/animal_metagenome/long-mg-dog/02_polishing/02_POLCA/D0${n}/D0${n}_PP1_PolcaCorr.fasta \
    -l ALL_SR,ALL_LR,ALL_multi -t 32 --write_bin_evals --write_bins \
    -o /data/anna/animal_metagenome/long-mg-dog/05_dereplication/00_dastool/D0${n}
  done
