#!/bin/bash

source ~/miniconda3/bin/activate root
conda activate hmmer2

BASE_DIR="/data/Projects/ShanghaiDogs"
RNAMMER_PATH="/data/anna/other_software/rnammer-1.2/rnammer"

DATA_DIR="$BASE_DIR/data/ShanghaiDogsMAGs"
OUTPUT_DIR="$BASE_DIR/intermediate-outputs/rnammer"
TMP_DIR="$BASE_DIR/intermediate-outputs/rnammer/tmp_rnammer"

mkdir -p "$TMP_DIR"  # tmp directory

cd "$DATA_DIR" || exit 1  # Exit if directory change fails

for file in *.fna.gz; do
  echo "Processing MAG: $file"
  temp_file="$TMP_DIR/$(basename "$file" .gz)"
  zcat "$file" > "$temp_file"
  cd "$TMP_DIR"
  $RNAMMER_PATH -S bac -m lsu,ssu,tsu -gff - "$temp_file" >> "$OUTPUT_DIR/rnammer_output.txt"
  rm "$temp_file"
  cd -
done
