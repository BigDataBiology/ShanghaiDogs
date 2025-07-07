# SHD SMORF Resource generation

## Step 1
```bash

Run_GMSC_Mapper.sh
```

Runs GMSC mapper on all samples in the SHD dataset.

## Step 2
```bash
python SHD_100AA_SmORFs.py
```

Removes redundancy and generates the 100AA version of the resource.

## Step 3

```bash
bash Cluster_SmORFs.sh
```

Runs cd-hit on the 100AA version to cluster the SmORFs at 90% AA identity.

## Step 4
```bash
python SHD_Clusters.py
```
Parses the output of cd-hit and generates the 90AA resource files


