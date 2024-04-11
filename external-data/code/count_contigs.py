import gzip
import os
import glob

# Define the directory where your gzipped files are located
genbank_dir = '/data/Projects/ShanghaiDogs/external-data/data/NCBI_genomes_ref/genbank/bacteria/'
refseq_dir = '/data/Projects/ShanghaiDogs/external-data/data/NCBI_genomes_ref/refseq/bacteria/'
output_path = '/data/Projects/ShanghaiDogs/external-data/data/NCBI_genomes_ref/contigs_count.txt'

# Initialize count variable
total_count = 0

# Loop through each folder in the base directory
with open(output_path,'w') as output_file:
    output_file.write("Filename,Number\n")
    for folder in os.listdir(genbank_dir):
        folder_path = os.path.join(genbank_dir, folder)
        if os.path.isdir(folder_path):
        # Loop through each gzipped file in the folder
            for file_path in glob.glob(os.path.join(folder_path, '*.gz')):
                count = 0
                print(file_path)
                with gzip.open(file_path, 'rt') as file:
                    for line in file:
                        if line.startswith('>'):
                            count += 1
                total_count += count
                # Write the filename and count to the output file
                filename = os.path.basename(file_path)
                output_file.write(f"{filename},{count}\n")
    for folder in os.listdir(refseq_dir):
        folder_path = os.path.join(refseq_dir, folder)
        if os.path.isdir(folder_path):
            # Loop through each gzipped file in the folder
            for file_path in glob.glob(os.path.join(folder_path, '*.gz')):
                count = 0
                print(file_path)
                with gzip.open(file_path, 'rt') as file:
                    for line in file:
                        if line.startswith('>'):
                            count += 1
                total_count += count
                # Write the filename and count to the output file
                filename = os.path.basename(file_path)
                output_file.write(f"{filename},{count}\n")


