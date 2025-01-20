library(dplyr)
library(ggplot2)
library(svglite)  # For saving plots as SVG

# Import maaslin2 results
# Reference categories were: "Study,This_study","env_classification,Dog Pet","Animal_age_simplified,Senior","Size_class,large","Sex,Male"
# Maaslin2 was run everytime with a cohort excluding unknowns/unclassified depending on which variable we were assessing 

maaslin_out <- read.delim("D:/AMAX_server/SingleM/maaslin_out_ENV_PetREF/significant_results.tsv", header = TRUE) # Pet as Ref
#maaslin_out <- read.delim("D:/AMAX_server/SingleM/maaslin_out_AGE_SeniorREF/significant_results.tsv", header = TRUE) # Senior as Ref
#maaslin_out <- read.delim("D:/AMAX_server/SingleM/maaslin_out_SIZE_LargeREF/significant_results.tsv", header = TRUE) # Large as Ref
#maaslin_out <- read.delim("D:/AMAX_server/SingleM/maaslin_out_SEX_MaleREF/significant_results.tsv", header = TRUE) # Male as Ref

category <- "Environment"
#category <- "Age"
#category <- "Size" # few significant if min qval 0.05
#category <- "Sex" # no significant

# Transform qvals for plotting
maaslin_out <- maaslin_out %>%
  mutate(`-log10(qvalue)` = -log10(qval))

# Create a species column
maaslin_out <- maaslin_out %>%
  mutate(species = sub(".*s__", "", feature))
maaslin_out <- maaslin_out %>%
  mutate(species = gsub("\\.", " ", species))

# Nicer metadata column nomenclature
maaslin_out <- maaslin_out %>%
  mutate(metadata = gsub("Animal_age_simplified", "Age", metadata),
         metadata = gsub("env_classification", "Environment", metadata),
         metadata = gsub("Size_class", "Size", metadata))

# Stricter filtering
maaslin_filt <- subset(maaslin_out,`-log10(qvalue)` >= -log10(0.01))
maaslin_filt <- subset(maaslin_filt,value != 'unknown')

# Add columns for absolute coefficient and sign of the coefficient
maaslin_filt <- maaslin_filt %>%
  mutate(abs_coef = abs(coef), 
         sign_coef = ifelse(coef > 0, "Positive", "Negative"))

# Identify the maximum `-log10(qvalue)` for each species and metadata feature
max_qvalue_per_feature <- maaslin_filt %>%
  group_by(species) %>%
  filter(`-log10(qvalue)` == max(`-log10(qvalue)`)) %>%
  ungroup()

# Add a column to indicate the metadata feature with the maximum `-log10(qvalue)`
maaslin_filt <- maaslin_filt %>%
  left_join(max_qvalue_per_feature %>% select(species, metadata) %>% rename(max_metadata = metadata), by = "species")

# Filter out species where the most significant feature is anything but the evaluated category
maaslin_filt <- maaslin_filt %>%
  filter(max_metadata == category)

# Value_order
value_order <- c("Other variables", "Wild Canid", "Dog Free_roaming", "Dog Colony")
#value_order <- c("Other variables", "Young", "Adult")
#value_order <- c("Other variables", "medium", "small")
#value_order <- c("Other variables", "female")

maaslin_filt <- maaslin_filt %>%
  mutate(value = factor(value, levels = value_order)) %>%
  arrange(value, `-log10(qvalue)`)

# Replace NA to other metadata variable
maaslin_filt <- maaslin_filt %>%
  mutate(value = as.character(value)) %>%  # Convert to character
  mutate(value = ifelse(is.na(value), "Other variables", value)) %>%  # Replace NA
  mutate(value = factor(value, levels = value_order))  # Convert back to factor

# Transform species to factor for the ordering
maaslin_filt <- maaslin_filt %>%
  mutate(species = factor(species, levels = unique(species)))

# Plot all together in a single plot - exclude species that study is the max

## ENVIRONMENT

p <- ggplot(maaslin_filt, aes(x = `-log10(qvalue)`, y = species, color = value, fill = value)) +
  # Circles
  geom_point(aes(size = abs_coef*2), shape = 21, alpha = 0.6, stroke = 1.2) +
  # Triangles inside circles
  geom_point(aes(size = 3, shape = sign_coef), fill = "white", color = "black") +
  geom_vline(xintercept = 2, linetype = "dashed", color = "darkgrey") +
  scale_shape_manual(values = c("Positive" = 24, "Negative" = 25)) + # triangles pointing up and down
  scale_color_manual(values = c("Wild Canid" = '#e7298a', "Dog Colony" = "#d95f02", 
                                "Dog Free_roaming" = "#e6ab02", "Other variables" = "lightgrey"),
                     breaks = rev(value_order)) +
  scale_fill_manual(values = c("Wild Canid" = '#e7298a', "Dog Colony" = "#d95f02", 
                               "Dog Free_roaming" = "#e6ab02", "Other variables" = "lightgrey"),
                    guide = "none") +
  scale_size_continuous(name = "Absolute Coef", 
                        guide = guide_legend(override.aes = list(shape = 21))) + 
  labs(title = "",  # Reference categories: Pet
       x = "-log10(qvalue)", y = "",
       color = "Metadata") +
  xlim(1, 8) + 
  theme_minimal() +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 1, size = 10), 
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10, face = "italic"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(color = "gray24", linewidth = 0.3))

#print(p)

# Save the plot as an SVG file
ggsave("D:/AMAX_server/SingleM/maaslin_figures/maaslin_ENV_pet.svg", plot = p,  width = 165, height = 60, units = "mm")

env_sign_sp <- unique(maaslin_filt$species[maaslin_filt$max_metadata == 'Environment'])
env_sign_sp <- as.character(env_sign_sp)
print(env_sign_sp)

## AGE
p <- ggplot(maaslin_filt, aes(x = `-log10(qvalue)`, y = species, color = value, fill = value)) +
  geom_point(aes(size = abs_coef, shape = sign_coef, stroke = 1.2), alpha = 0.6) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "darkgrey") +
  scale_shape_manual(values = c("Positive" = 24, "Negative" = 25)) +  # Filled circles for positive, empty circles for negative
  scale_color_manual(values = c("Young" = '#e6ab02', "Adult" = "#d95f02", "Other variables" = "lightgrey"),
                     breaks = rev(value_order)) +
  scale_fill_manual(values = c("Young" = '#e6ab02', "Adult" = "#d95f02", "Other variables" = "lightgrey"),
                    guide = "none") +
  scale_size_area(name = "Absolute Coef", 
                  max_size = 6,  # Adjust max size for better visual contrast
                  guide = guide_legend(override.aes = list(shape = 24))) + 
  labs(title = "", 
       x = "-log10(qvalue)",
       y = "",
       color = "Metadata",
       size = "Absolute Coef") +
  xlim(1.8, 5.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(face = "italic"),
        text = element_text(size = 10))  # Set all font sizes to 10
#print(p)

# Save the plot with the specified size in mm
ggsave("D:/AMAX_server/SingleM/maaslin_figures/maaslin_AGE_senior.svg", plot = p, width = 130, height = 100, units = "mm", device = "svg")

age_sign_sp <- unique(maaslin_filt$species[maaslin_filt$max_metadata == 'Age'])
age_sign_sp <- as.character(age_sign_sp)
print(age_sign_sp)

## AGE with YOUNG as a reference
maaslin_out <- read.delim("D:/AMAX_server/SingleM/maaslin_out_AGE_Young/significant_results.tsv", header = TRUE) # Young as Ref
category <- "Age"

# Transform qvals and create species column
maaslin_out <- maaslin_out %>%
  mutate(`-log10(qvalue)` = -log10(qval),
         species = sub(".*s__", "", feature) %>% gsub("\\.", " ", .))

# Nicer metadata column nomenclature
maaslin_out <- maaslin_out %>%
  mutate(metadata = gsub("Animal_age_simplified", "Age", metadata),
         metadata = gsub("env_classification", "Environment", metadata),
         metadata = gsub("Size_class", "Size", metadata))

# Filter out low q-values and unknowns
maaslin_filt <- subset(maaslin_out,`-log10(qvalue)` >= -log10(0.01))
maaslin_filt <- subset(maaslin_filt,value != 'unknown')

# Add columns for absolute coefficient and sign of the coefficient
maaslin_filt <- maaslin_filt %>%
  mutate(abs_coef = abs(coef), 
         sign_coef = ifelse(coef > 0, "Positive", "Negative"))

# Identify the most significant metadata feature per species
max_qvalue_per_feature <- maaslin_filt %>%
  group_by(species) %>%
  filter(`-log10(qvalue)` == max(`-log10(qvalue)`)) %>%
  ungroup()

maaslin_filt <- maaslin_filt %>%
  left_join(max_qvalue_per_feature %>% select(species, metadata) %>% rename(max_metadata = metadata), by = "species") %>%
  filter(max_metadata == "Age")

# Ensure 'value' column is properly set as a factor with correct levels
value_order <- c("Senior","Adult", "Other variables")
maaslin_filt <- maaslin_filt %>%
  mutate(value = factor(value, levels = value_order))  # Convert 'value' to factor

# Define species sorting and plot
maaslin_filt <- maaslin_filt %>%
  arrange(desc(`-log10(qvalue)`)) %>%
  mutate(species = factor(species, levels = rev(unique(species))))

# Replace NA to other metadata variable
maaslin_filt <- maaslin_filt %>%
  mutate(value = as.character(value)) %>%  # Convert to character
  mutate(value = ifelse(is.na(value), "Other variables", value)) %>%  # Replace NA
  mutate(value = factor(value, levels = value_order))  # Convert back to factor

# Plot
p <- ggplot(maaslin_filt, aes(x = `-log10(qvalue)`, y = species, color = value, fill = value)) +
  geom_point(aes(size = abs_coef, shape = sign_coef, stroke = 1.2), alpha = 0.6) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "darkgrey") +
  scale_shape_manual(values = c("Positive" = 24, "Negative" = 25)) +  # Filled circles for positive, empty circles for negative
  scale_color_manual(values = c("Senior" = '#e6ab02', "Adult" = "#d95f02", "Other variables" = "lightgrey"),
                     breaks = rev(value_order)) +
  scale_fill_manual(values = c("Senior" = '#e6ab02', "Adult" = "#d95f02", "Other variables" = "lightgrey"),
                    guide = "none") +
  scale_size_area(name = "Absolute Coef", 
                  max_size = 6,  # Adjust max size for better visual contrast
                  guide = guide_legend(override.aes = list(shape = 24))) + 
  labs(x = "-log10(qvalue)", y = "", color = "Metadata", size = "Absolute Coef") +
  xlim(1.8, 5.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(face = "italic"),
        text = element_text(size = 10))  # Set all font sizes to 10

#print(p)

# Save the plot with the specified size in mm
ggsave("D:/AMAX_server/SingleM/maaslin_figures/maaslin_AGE_young.svg", plot = p, width = 130, height = 100, units = "mm", device = "svg")

