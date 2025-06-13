import os
from pathlib import Path

# paths
base_dir = "/work/microbiome/shanghai_dogs/intermediate-outputs/GMSC_MAPPER"
output_file = os.path.join(base_dir, "All_mapped.smorfs.faa")

sample_dirs = [
    "D000_PP1_PolcaCorr", "D001_PP1_PolcaCorr", "D002_PP1_PolcaCorr", "D003_PP1_PolcaCorr",
    "D004_PP1_PolcaCorr", "D005_PP1_PolcaCorr", "D006_PP1_PolcaCorr", "D007_PP1_PolcaCorr",
    "D008_PP1_PolcaCorr", "D010_PP1_PolcaCorr", "D011_PP1_PolcaCorr", "D012_PP1_PolcaCorr",
    "D013_PP1_PolcaCorr", "D014_PP1_PolcaCorr", "D015_PP1_PolcaCorr", "D016_PP1_PolcaCorr",
    "D017_PP1_PolcaCorr", "D018_PP1_PolcaCorr", "D019_PP1_PolcaCorr", "D020_PP1_PolcaCorr",
    "D021_PP1_PolcaCorr", "D022_PP1_PolcaCorr", "D023_PP1_PolcaCorr", "D024_PP1_PolcaCorr",
    "D025_PP1_PolcaCorr", "D026_PP1_PolcaCorr", "D027_PP1_PolcaCorr", "D028_PP1_PolcaCorr",
    "D029_PP1_PolcaCorr", "D030_PP1_PolcaCorr", "D031_PP1_PolcaCorr", "D032_PP1_PolcaCorr",
    "D033_PP1_PolcaCorr", "D034_PP1_PolcaCorr", "D035_PP1_PolcaCorr", "D036_PP1_PolcaCorr",
    "D037_PP1_PolcaCorr", "D038_PP1_PolcaCorr", "D039_PP1_PolcaCorr", "D040_PP1_PolcaCorr",
    "D041_PP1_PolcaCorr", "D042_PP1_PolcaCorr", "D043_PP1_PolcaCorr", "D044_PP1_PolcaCorr",
    "D045_PP1_PolcaCorr", "D046_PP1_PolcaCorr", "D047_PP1_PolcaCorr", "D048_PP1_PolcaCorr",
    "D049_PP1_PolcaCorr", "D050_PP1_PolcaCorr", "D051_PP1_PolcaCorr", "D052_PP1_PolcaCorr"
]

# Counter for smORF IDs
smorf_count = 0

with open(output_file, "w") as outfile:
    # Iterate through each sample directory
    for sample in sample_dirs:
        sample_path = os.path.join(base_dir, sample, "mapped.smorfs.faa")
        
        # process FASTA file
        with open(sample_path, "r") as infile:
            sequence = ""
            for line in infile:
                line = line.strip()
                if line.startswith(">"):
                    # Write previous sequence if exists
                    if sequence:
                        # Format ID as SHD_SM.XXX.XXX
                        first_part = smorf_count // 1000
                        second_part = smorf_count % 1000
                        new_id = f"SHD_SM.{str(first_part).zfill(3)}.{str(second_part).zfill(3)}"
                        outfile.write(f">{new_id}\n{sequence}\n")
                        smorf_count += 1
                        sequence = ""
                else:
                    # Accumulate sequence lines
                    sequence += line
            #last sequence
            if sequence:
                first_part = smorf_count // 1000
                second_part = smorf_count % 1000
                new_id = f"SHD_SM.{str(first_part).zfill(3)}.{str(second_part).zfill(3)}"
                outfile.write(f">{new_id}\n{sequence}\n")
                smorf_count += 1

print(f"Generated {output_file} with {smorf_count} smORFs.")
