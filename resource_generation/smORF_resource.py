conda activate gmscmapper

# paths
INPUT_DIR='shanghai_dogs\data\ShanghaiDogsAssemblies'
OUTPUT_DIR='shanghai_dogs/intermediate-outputs/GMSC_MAPPER_TEST'
GMSC='shanghai_dogs/intermediate-outputs/GMSC_Mapper_dogs/GMSC-mapper'

# Loop over all files with .fa, .fasta, or .fna extensions in the input directory
for file in "$INPUT_DIR"/*.{fa,fasta,fna}; do
    # Skip if no files are found
    [ -e "$file" ] || continue

    filename=$(basename "$file")
    filename_no_ext="${filename%.*}"
    outdir="$OUTPUT_DIR/$filename_no_ext"

    if [ -d "$outdir" ] && [ "$(ls -A "$outdir")" ]; then
        echo "✅ Skipping $filename (output exists)"
        continue
    fi

    echo "🚀 Running gmsc-mapper on $filename"
    mkdir -p "$outdir"
    gmsc-mapper -i "$file" -o "$outdir" --dbdir "$GMSC/db" -t 24
done

