import pandas as pd


qual_report = pd.read_csv('data/ShanghaiDogsTables/SHD_bins_MIMAG_report.csv')
NCBI_ref = qual_report['GTDBtk fastani Ref'].dropna().unique().tolist()
print(NCBI_ref)

out_path = "external-data/data/GTDBtk_fastani_ref_genomes/fastani_ref_genome_assemblies.txt"

with open(out_path, 'w') as file:
    for assembly in NCBI_ref:
        file.write(str(assembly) + '\n')
