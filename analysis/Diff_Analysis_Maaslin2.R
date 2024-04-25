library(Maaslin2)
library(dplyr)
library(tidyverse)
library(phyloseq)
library(vegan)

setwd("C:/Users/Anna/BDB-Anna-work/8-Pet_microbiome_long/04_R_ANALYSIS/analysis")

# Import data
OTU_tab <- read.csv("REP_OTU_count_filt_S3.5.rib_prot_S2_rpsB.csv", header=TRUE, sep = ",", row.names=1)
OTU_tab <- t(OTU_tab)
metadata <- read.csv("metadata_filt.csv", header = TRUE, sep = ",",row.names = 1)

# Check categories & create a strict subsets

metadata %>% dplyr::count(Study)
metadata_filt <- metadata %>% filter(!str_detect(metadata$Study, "Rampelli_2021_coprolite"))
metadata_filt <- metadata_filt %>% filter(!str_detect(metadata_filt$Study, "Pehrsson_2016"))
metadata_filt <- metadata_filt %>% filter(!str_detect(metadata_filt$Study, "Alessandri_2019_dogs"))

metadata_filt %>% dplyr::count(env_classification)
metadata_filt <- metadata_filt %>% filter(!str_detect(metadata_filt$env_classification, "Dog others"))
metadata_filt <- metadata_filt %>% filter(!str_detect(metadata_filt$env_classification, "Dog Ancient"))

OTU_filt <- OTU_tab[row.names(OTU_tab) %in% row.names(metadata_filt) ,]
OTU_filt <- OTU_filt[, which(colSums(OTU_filt) != 0)] #filter motus that have 0 counts

### Run Maaslin2 in OTU tables

fit_data = Maaslin2(
  input_data = OTU_filt,
  input_metadata = metadata_filt,
  min_prevalence = 0.1,
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  output = "maaslin_out/env_classification_study",
  fixed_effects = c("env_classification",'Study'),
  reference = c("env_classification,Dog Pet","Study,This_study"))
