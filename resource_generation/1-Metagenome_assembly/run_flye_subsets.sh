#!/bin/bash

source ~/miniconda3/bin/activate root
conda activate Flye-2.9.2

ont_subsets="/data/anna/animal_metagenome/long-mg-dog/00_quality_control/01_trim_filter_data/ont/porechop/subsets"
out_flye="/data/anna/animal_metagenome/long-mg-dog/01_assembly/subset_assemblies"

# Flye assembly

# Samples present different raw error rates - adjusted accordingly
declare -a list_03=("07" "10" "11" "12" "13" "14" "15" "16" "17" "19" "21" "22" "25" "27" "29" "30" "32" "34" "35" "38" "39" "40" "42" "44" "45" "46" "47" "48" "50" "51")
declare -a list_04=("00" "01" "03" "06" "08" "18" "20" "23" "26" "31" "36" "37" "41" "43" "49" "52" )
declare -a list_05=("02" "04" "05" "24" "28" "33" "46")

for n in ${list_03[@]}
do
    echo "Processing D0${n}..."
    for subset in 10G 20G
    do
        input="${ont_subsets}/D0${n}_${subset}_filt_porechop.fastq.gz"
        output="${out_flye}/D0${n}_${subset}/"

        if [[ -f "$input" ]]; then
            flye -o "$output" --threads 32 --nano-hq "$input" --meta --read-error 0.05
        else
            echo "File $input not found, skipping..."
        fi
    done
done

for n in ${list_04[@]}
do
    echo "Processing D0${n}..."
    for subset in 10G 20G
    do
        input="${ont_subsets}/D0${n}_${subset}_filt_porechop.fastq.gz"
        output="${out_flye}/D0${n}_${subset}/"

        if [[ -f "$input" ]]; then
            flye -o "$output" --threads 32 --nano-hq "$input" --meta --read-error 0.04
        else
            echo "File $input not found, skipping..."
        fi
    done
done


for n in ${list_05[@]}
do
    echo "Processing D0${n}..."
    for subset in 10G 20G
    do
        input="${ont_subsets}/D0${n}_${subset}_filt_porechop.fastq.gz"
        output="${out_flye}/D0${n}_${subset}/"

        if [[ -f "$input" ]]; then
            flye -o "$output" --threads 32 --nano-hq "$input" --meta --read-error 0.05
        else
            echo "File $input not found, skipping..."
        fi
    done
done
