import os
import csv
from collections import defaultdict
import fasta

# Paths
base_dir = "/work/microbiome/shanghai_dogs/intermediate-outputs/GMSC_MAPPER"
input_fasta = os.path.join(base_dir, "100AA_SmORFs_sequences.faa")
cluster_tsv = os.path.join(base_dir, "SHD_Clusters.tsv")
metadata_100aa_tsv = os.path.join(base_dir, "100AA_SmORFs_metadata.tsv")
habitat_taxonomy_100aa_tsv = os.path.join(base_dir, "100AA_SmORFs_habitat_taxonomy.tsv")
output_fasta = os.path.join(base_dir, "90AA_SmORFs_sequences.faa")
output_metadata_tsv = os.path.join(base_dir, "90AA_SmORFs_metadata.tsv")
output_habitat_taxonomy_tsv = os.path.join(base_dir, "90AA_SmORFs_habitat_taxonomy.tsv")


