import os
import re

# Base directory for GMSC-Mapper output
GMSC_output_data = "shanghai_dogs/intermediate-outputs/GMSC_MAPPER_TEST"
sample_dirs = sorted([d for d in os.listdir(GMSC_output_data) if re.match(r"D\d{3}_PP1_PolcaCorr", d)])

for sample_dir in sample_dirs:
    # Input and output FASTA paths
    input_faa = os.path.join(GMSC_output_data, sample_dir, "mapped.smorfs.faa")
    output_faa = os.path.join(GMSC_output_data, sample_dir, "SHD_mapped_smorfs.faa")

    # Extract sample number from directory name
    sample_number = re.match(r"D(\d+)_PP1_PolcaCorr", sample_dir).group(1)

    # Process FASTA: rename headers to SHD_SM.<sample>.<count>
    with open(input_faa, "r") as f_in, open(output_faa, "w") as f_out:
        smorf_count = 1
        for line in f_in:
            if line.startswith(">"):
                f_out.write(f">SHD_SM.{sample_number}.{smorf_count:05d}\n")
                smorf_count += 1
            else:
                f_out.write(line)

    print(f"[DONE] {sample_dir}: Renamed {smorf_count - 1} smORFs → {output_faa}")
