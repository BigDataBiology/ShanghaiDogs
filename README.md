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

1. `ShanghaiDogsFastQ/`: Raw FastQs from service provider. FQs will be uploaded, Fast5 maybe
2. `ShanghaiDogsMetadata/`: Metadata for the samples
3. `ShanghaiDogsAssemblies/`: Polished assemblies
4. `ShanghaiDogsMAGs/`: MAGs (FASTA files)
5. `ShanghaiDogsMAGAnnotations/`: Annotations of MAGs, including `EMapper` and `Barrnap` subdirectories
6. `ShanghaiDogsMAGsTables/`: Tables with MAGs, including `singleM` and `mOTU` subdirectories
7. `ShanghaiDogs_OtherResources/`: Other resources, including gene and smORF catalogues (and their respective annotations)

### External data

1. GMGC MAGs
2. singleM profiles from sandpiper (inside `dog_microbiome_archive_otu_tables` subdirectory)
3. mOTU profiles for canid/felid

### Preprocessed

1. Quality controled FQs
2. Not fully polished assemblies
3. singleM/mOTU output tables

