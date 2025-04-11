#!/bin/bash

source ~/miniconda3/bin/activate root
conda activate R-env

# Main paths
ont_data_1='/data/anna/animal_metagenome/long-mg-dog/00_quality_control/00_raw_clean_data/ont/01_seq_summaries'
ont_data_2='/data/Public/Result-X101SC23043778-Z01-J005-B1-6_delivery_20230818/ONT'
output='/data/anna/animal_metagenome/long-mg-dog/00_quality_control/00_raw_clean_data/ont/02_MinIONQC_OK'

# FIRST BATCH OF SAMPLES (D000-D008)
# Moved and renamed the sequencing summary files, D00X.sequencing_summary.txt not read by MinIONQC.R
for n in {00..08}
  do
    Rscript ~/MinIONQC.R -i ${ont_data_1}/D0${n}/sequencing_summary.txt -q 10 -o ${output}
  done

# D001
cd ~/animal_metagenome/long-mg-dog/00_quality_control/00_raw_clean_data/ont/01_seq_summaries/D001
head -n 1 sequencing_summary_0.txt > D001_head.txt
sed -i '1d' sequencing_summary_0.txt
sed -i '1d' sequencing_summary_1.txt
sed -i '1d' sequencing_summary_2.txt
cat D001_head.txt sequencing_summary_0.txt sequencing_summary_1.txt sequencing_summary_2.txt > sequencing_summary.txt
Rscript ~/MinIONQC.R -i sequencing_summary.txt -q 10 -o ${output}/D001 &
rm sequencing_summary_0.txt sequencing_summary_1.txt sequencing_summary_2.txt D001_head.txt

# D006
cd ~/animal_metagenome/long-mg-dog/00_quality_control/00_raw_clean_data/ont/01_seq_summaries/D006
head -n 1 sequencing_summary_0.txt > D006_head.txt
sed -i '1d' sequencing_summary_0.txt
sed -i '1d' sequencing_summary_1.txt
sed -i '1d' sequencing_summary_2.txt
cat D006_head.txt sequencing_summary_0.txt sequencing_summary_1.txt sequencing_summary_2.txt > sequencing_summary.txt
Rscript ~/MinIONQC.R -i sequencing_summary.txt -q 10 -o ${output}/D006 &
rm sequencing_summary_0.txt sequencing_summary_1.txt sequencing_summary_2.txt D006_head.txt

# NEW SAMPLES (D010-D052)
for n in {10..52}
  do
    Rscript ~/MinIONQC.R -i ${ont_data_2}/D0${n}/sequencing_summary* -q 10 -o ${output}
  done

# Create QC_all report
cd ${output}

for n in {00..52}; do
  if [ -d "D0${n}" ]; then
    echo "D0${n}" >> sample_id
    sed -n '30p' D0${n}/summary.yaml >> Gbases
    sed -n '32p' D0${n}/summary.yaml >> N50
    sed -n '34p' D0${n}/summary.yaml >> median_length
    sed -n '37p' D0${n}/summary.yaml >> median_q
  else
    echo "D0${n} directory does not exist."
  fi
done

paste sample_id Gbases N50 median_length median_q > /data/Projects/ShanghaiDogs/intermediate-outputs/00_quality_control/QC_ONT_raw_reads.txt

# Clean up temporary files
rm sample_id Gbases N50 median_q median_length