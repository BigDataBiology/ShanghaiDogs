# SHD SMORF Resource Generation

This pipeline generates small Open Reading Frame (SmORF) resources from the SHD dataset

### Step 1: GMSC Mapping
```bash
Run_GMSC_Mapper.sh
```
Executes GMSC mapper on all samples in the SHD dataset. The GMSC mapper identifies the small open reading frames across all the samples.

**Tool:** [GMSC-mapper](https://github.com/BigDataBiology/GMSC-mapper) with default usage parameters

### Step 2: Redundancy Removal
```bash
python SHD_100AA_SmORFs.py
```
Processes the GMSC mapper output to remove redundant sequences and generates the 100 amino acid (AA) version of the SmORF resource. This step creates a non-redundant set of SmORFs for downstream analysis. This creates files two file namely 100AA_SmORFs_origins.tsv.gz ; 100AA_SmORFs_sequences.faa.gz 

### Step 3: Sequence Clustering
```bash
bash Cluster_SmORFs.sh
```
Performs sequence clustering using CD-HIT on the 100AA SmORF sequences. Clusters are formed based on 90% amino acid identity using default CD-HIT parameters.

**Tool:** [CD-HIT](https://github.com/weizhongli/cdhit) with default parameters

### Step 4: Cluster Mapping
```bash
python SHD_Clusters.py
```
Parses the CD-HIT clustering output and maps 100AA SmORFs to 90AA cluster representatives.

**PARAMETERS**

- **SmORF ID**: Unique identifier for each small Open Reading Frame.
- **Sample ID**: Identifier of the sample from which the SmORF was derived.
- **Contig**: Contig name where the SmORF is located.
- **Coordinates**: Genomic coordinates of the SmORF on the contig.
- **Strand**: DNA strand orientation (forward/reverse) of the SmORF.

