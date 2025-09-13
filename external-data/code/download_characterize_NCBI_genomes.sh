#!/bin/bash


BASEDIR=$PWD
cd ${BASEDIR}/external-data/data/NCBI_genomes_ref/

ncbi-genome-download --assembly-accessions fastani_ref_genome_assemblies.txt --formats 'fasta' 'bacteria' #GCF assemblies
ncbi-genome-download --assembly-accessions fastani_ref_genome_assemblies.txt --formats 'fasta' 'bacteria' -s 'genbank' #GCA_ assemblies

# Count contigs of the genomes

cd ${BASEDIR}/external-data/code
python count_contigs.py

# Run checkm2 into the genomes

conda deactivate
conda activate checkm2

cd ${BASEDIR}/external-data/data/NCBI_genomes_ref/
checkm2 predict -x fna.gz --threads 32 --input */bacteria/*/GC* \
--output-directory checkm2/

# Run barrnap into the genomes

conda deactivate
conda activate barrnap

cd ${BASEDIR}/external-data/data/NCBI_genomes_ref/
gzip -d */bacteria/*/GC*.fna.gz

path_genbank=${BASEDIR}/external-data/data/NCBI_genomes_ref/genbank/bacteria/
path_refseq=${BASEDIR}/external-data/data/NCBI_genomes_ref/refseq/bacteria/

for file in "$path_genbank"*/GC*.fna
  do
    filename=$(basename "$file" _genomic.fna)
    GC_name=$(echo "$filename" | cut -d'_' -f1,2)
    echo ${GC_name}
    barrnap --kingdom 'bac' --threads 24 ${file} \
    --outseq barrnap/fasta/${GC_name}_ribosomal.fa > \
    barrnap/out/${GC_name}_barrnap.txt
  done

for file in "$path_refseq"*/GC*.fna
  do
    filename=$(basename "$file" _genomic.fna)
    GC_name=$(echo "$filename" | cut -d'_' -f1,2)
    echo ${GC_name}
    barrnap --kingdom 'bac' --threads 24 ${file} \
    --outseq barrnap/fasta/${GC_name}_ribosomal.fa > \
    barrnap/out/${GC_name}_barrnap.txt
  done

# Run tRNA-scan into the genomes

conda deactivate
conda activate trnascan-se

for file in "$path_genbank"*/GC*.fna
  do
    filename=$(basename "$file" _genomic.fna)
    GC_name=$(echo "$filename" | cut -d'_' -f1,2)
    echo ${GC_name}
    tRNAscan-SE -B -o tRNAs/${GC_name}_trna.out ${file}
  done

for file in "$path_refseq"*/GC*.fna
  do
    filename=$(basename "$file" _genomic.fna)
    GC_name=$(echo "$filename" | cut -d'_' -f1,2)
    echo ${GC_name}
    tRNAscan-SE -B -o tRNAs/${GC_name}_trna.out ${file}
  done

gzip */bacteria/*/GC*.fna
