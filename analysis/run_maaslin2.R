library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(DESeq2)
library(tidyverse)
library(Maaslin2)
library(dplyr)

#Import metadata

metadata <- read.delim("ShanghaiDogs/data/ShanghaiDogsMetadata/REP_canid_metadata.csv", header = TRUE, sep = ',')
rownames(metadata) <- metadata[, 1]
metadata <- metadata[, -1]
metadata$env_classification <- gsub(" captive","",metadata$env_classification)
metadata$env_classification <- gsub("Dog others","Dog undet",metadata$env_classification)

# Import otu/taxonomy table
otutab <- read.delim("ShanghaiDogs/intermediate-outputs/singlem-profiling/tax-profiles/all-dog-tax-profile-species.tsv", header = TRUE, sep = '\t')

# Import list of samples to be included for the analysis
# Are the samples considered for beta diversity analysis
ls_samples <- scan("D:/AMAX_server/SingleM/samples_ls.csv", sep=",", what = "")
ls_samples <- as.list(ls_samples)

# create taxonomy table
taxtab <- as.data.frame(otutab$taxonomy)
colnames(taxtab) <- "taxonomy"
rownames(taxtab) <- taxtab$taxonomy
taxtab <- separate(taxtab, taxonomy, into = c("Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "; ")
taxtab[is.na(taxtab)] <- "unassigned"
taxtab <- taxtab[, -1] # remove first column
taxtab_matrix <- as.matrix(taxtab)
# write.csv(taxtab,"ShanghaiDogs/intermediate-outputs/singlem-profiling/Taxonomy-all-dog-tax.csv")

# Pre-process otu/taxonomy-profile: update naming, remove underscores, etc.
SRR_archived <- read.delim("ShanghaiDogs/external-data/data/dog_microbiome_archive_otu_tables/SRR_Acc_List_Metadata.txt", header = TRUE, sep = '|')
SRR_archived <- SRR_archived[-1, ] # remove first row
SRR_archived$run <- gsub(" ", "", SRR_archived$run)
SRR_archived$biosample <- gsub(" ", "", SRR_archived$biosample)

for (i in seq_along(colnames(otutab))) {
  if (colnames(otutab)[i] %in% SRR_archived$run) {
    colnames(otutab)[i] <- SRR_archived$biosample[SRR_archived$run == colnames(otutab)[i]]
  }
}

rownames(otutab) <- otutab[, 1] 
otutab <- otutab[, -1] # remove first column
colnames(otutab) <- gsub("_.*", "",colnames(otutab))

# Filter metadata and samples to match sample_ls and variable: ENV_CLASSIFICATION
metadata_filt <- metadata[rownames(metadata) %in% ls_samples,]
metadata_filt <- metadata_filt %>% filter(!str_detect(metadata_filt$env_classification, "Dog undet"))
metadata_filt <- metadata_filt %>% filter(!str_detect(metadata_filt$env_classification, "Dog Shelter"))

metadata_filt %>% dplyr::count(Study)
metadata_filt %>% dplyr::count(env_classification)

OTU_filt <- otutab[,colnames(otutab) %in% row.names(metadata_filt)]
OTU_filt <- OTU_filt[which(rowSums(OTU_filt) != 0), ] #filter otus that have 0 counts

fit_data = Maaslin2(
  input_data = OTU_filt, 
  input_metadata = metadata_filt, 
  min_prevalence = 0.1,
  normalization = "NONE", #already RA
  transform = "LOG",
  analysis_method = "LM",
  output = "ShanghaiDogs/intermediate-outputs/singlem-profiling/Maaslin2_output/maaslin_out_ENV_PetREF/", 
  fixed_effects = c("Study","env_classification","Size_class","Animal_age_simplified","Sex"),
  reference = c("Study,This_study","env_classification,Dog Pet","Size_class,large","Animal_age_simplified,Senior","Sex,Male"))

# Filter metadata and samples to match sample_ls and variable: AGE
metadata_filt <- metadata[rownames(metadata) %in% ls_samples,]
metadata_filt <- metadata_filt %>% filter(!str_detect(metadata_filt$env_classification, "Wild Canid"))
metadata_filt <- metadata_filt %>% filter(!str_detect(metadata_filt$Animal_age_simplified, "unknown"))
metadata_filt %>% dplyr::count(Animal_age_simplified)

OTU_filt <- otutab[,colnames(otutab) %in% row.names(metadata_filt)]
OTU_filt <- OTU_filt[which(rowSums(OTU_filt) != 0), ] #filter otus that have 0 counts

fit_data = Maaslin2(
  input_data = OTU_filt, 
  input_metadata = metadata_filt, 
  min_prevalence = 0.1,
  normalization = "NONE", #already RA
  transform = "LOG",
  analysis_method = "LM",
  output = "ShanghaiDogs/intermediate-outputs/singlem-profiling/Maaslin2_output/maaslin_out_AGE_seniorREF/", 
  fixed_effects = c("Study","env_classification","Size_class","Animal_age_simplified","Sex"),
  reference = c("Study,This_study","env_classification,Dog Pet","Size_class,large","Animal_age_simplified,Senior","Sex,Male"))

# Filter metadata and samples to match sample_ls and variable: SIZE
metadata_filt <- metadata[rownames(metadata) %in% ls_samples,]
metadata_filt <- metadata_filt %>% filter(!str_detect(metadata_filt$env_classification, "Wild Canid"))
metadata_filt <- metadata_filt %>% filter(!str_detect(metadata_filt$Size_class, "unclass"))
metadata_filt %>% dplyr::count(Size_class)

OTU_filt <- otutab[,colnames(otutab) %in% row.names(metadata_filt)]
OTU_filt <- OTU_filt[which(rowSums(OTU_filt) != 0), ] #filter otus that have 0 counts

fit_data = Maaslin2(
  input_data = OTU_filt, 
  input_metadata = metadata_filt, 
  min_prevalence = 0.1,
  normalization = "NONE", #already RA
  transform = "LOG",
  analysis_method = "LM",
  output = "ShanghaiDogs/intermediate-outputs/singlem-profiling/Maaslin2_output/maaslin_out_SIZE_LargeREF/", 
  fixed_effects = c("Study","env_classification","Size_class","Animal_age_simplified","Sex"),
  reference = c("Study,This_study","env_classification,Dog Pet","Size_class,large","Animal_age_simplified,Senior","Sex,Male"))


# Filter metadata and samples to match sample_ls and variable: SEX
metadata_filt <- metadata[rownames(metadata) %in% ls_samples,]
metadata_filt <- metadata_filt %>% filter(!str_detect(metadata_filt$env_classification, "Wild Canid"))
metadata_filt <- metadata_filt %>% filter(!str_detect(metadata_filt$Sex, "Unknown"))
metadata_filt <- metadata_filt %>% filter(!str_detect(metadata_filt$Sex, "unknown"))
metadata_filt %>% dplyr::count(Sex)

OTU_filt <- otutab[,colnames(otutab) %in% row.names(metadata_filt)]
OTU_filt <- OTU_filt[which(rowSums(OTU_filt) != 0), ] #filter otus that have 0 counts

fit_data = Maaslin2(
  input_data = OTU_filt, 
  input_metadata = metadata_filt, 
  min_prevalence = 0.1,
  normalization = "NONE", #already RA
  transform = "LOG",
  analysis_method = "LM",
  output = "ShanghaiDogs/intermediate-outputs/singlem-profiling/Maaslin2_output/maaslin_out_SEX_MaleREF/", 
  fixed_effects = c("Study","env_classification","Size_class","Animal_age_simplified","Sex"),
  reference = c("Study,This_study","env_classification,Dog Pet","Size_class,large","Animal_age_simplified,Senior","Sex,Male"))


