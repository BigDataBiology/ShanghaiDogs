# Shanghai Dogs


## Directories

- `data/` (not on git): raw data and results to upload (to ENA/Zenodo/...)
- `external-data/code` (on git): code to download external data
- `external-data/data` (not on git): external data downloaded
- `intermediate-outputs/` (not on git): results of preprocessing for convenience
- `resource_generation/`: scripts to generate assemblies, MAGs, annotations, ...
- `analysis/`
- `notes/`


### Data

1. Raw FastQs from service provider. FQs will be uploaded, Fast5 maybe
2. Metadata
3. Polished assemblies
4. MAGs
5. Annotation tables

### External data

1. GMGC MAGs
2. singleM profiles from sandpiper (inside `dog_microbiome_archive_otu_tables` subdirectory)
3. mOTU profiles for canid/felid

### Preprocessed

1. Quality controled FQs
2. Not fully polished assemblies
3. singleM/mOTU output tables

