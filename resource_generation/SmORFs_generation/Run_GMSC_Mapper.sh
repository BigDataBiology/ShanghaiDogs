#Install GMSC-mapper: https://github.com/BigDataBiology/GMSC-mapper
#!/bin/bash -l
#PBS -N gmsc_mapper_job
#PBS -l select=1:ncpus=24:mem=80gb
#PBS -l walltime=48:00:00
#PBS -m abe
#PBS -q microbiome

cd "$PBS_O_WORKDIR"

source ~/.bashrc
conda activate gmscmapper
#paths
input="/shanghai_dogs/data/ShanghaiDogsAssemblies"
output="/shanghai_dogs/intermediate-outputs/GMSC_MAPPER"
GMSC="/home/n12228516/tools/GMSC-mapper_git_version"

logfile="gmsc_errors.log"
summaryfile="gmsc_summary_report.log"

echo "GMSC Mapper error log - $(date)" > "$logfile"
echo "GMSC Mapper summary report - $(date)" > "$summaryfile"

find "$input" -maxdepth 1 -type f \( -name "*.fna" -o -name "*.fna.gz" \) | while read -r file; do
    if [ ! -e "$file" ]; then
        echo " File not found: $file at $(date)" >> "$logfile"
        continue
    fi

    filename=$(basename "$file")
    filename_no_ext="${filename%.fna}"
    filename_no_ext="${filename_no_ext%.gz}"

    outdir="$output/$filename_no_ext"
    mkdir -p "$outdir"

    echo "🚀 Running gmsc-mapper on $filename"
    if ! gmsc-mapper -i "$file" -o "$outdir" --dbdir "$GMSC/db" -t 24; then
        echo " GMSC-mapper execution failed for $filename at $(date)" >> "$logfile"
        echo "$filename: FAILED during execution" >> "$summaryfile"
        continue
    fi

    echo "$filename:Processed successfully" >> "$summaryfile"
done
