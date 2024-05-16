# This script runs repeated measures ANOVA on normalized (Normal - Scramble) PLI values for a 2x3 (Species: Human/Monkey X Category: Body/Face/Object) design 
# clean environment 
rm(list = ls())

# required packages 
library(readxl)
library(tidyverse)
library(car)
library(broom)
library(knitr)
library(rstatix)
library(dplyr)
library(ggpubr)
library(xlsx)
library(ggplot2)

# Define parent directory
dir_parent = "/Users/jane_chesley/Documents/Github/Dataset_dynamic/EEG-analysis/statistics/R"

# data to analyze 

# List all folders in the parent directory
# Do not list subfolders ('recursive = FALSE')
all_folders <- list.dirs(dir_parent, full.names = FALSE, recursive = FALSE)

# Filter folder names: keep only folders that have "PLI_normalized" in their name
filtered_folders <- all_folders[grep("PLI_normalized", all_folders, fixed = TRUE)]

# Construct file names from folder names 
file_names <- paste(filtered_folders, ".xlsx", sep = "")

# data to analyze 
all_data <- data.frame(folder_name = filtered_folders, file_name = file_names)

# run analysis for each data file 
for (i in 1:nrow(all_data)) {
  current_folder <- all_data[i, "folder_name"]
  current_file <- all_data[i, "file_name"]
  
  
  # Set working directory 
  setwd(paste(dir_parent, current_folder, sep = "/"))
  
  # Get input data in long format
  long_data <- read_excel(current_file)
  
  
  
  
  # Convert numerical to character labels 
  
  # Variable 1
  species_labels <- c("human", "monkey")
  # Replace numerical values with character labels
  long_data$Species <- species_labels[long_data$Species]
  
  # Variable 2
  category_labels <- c("body", "face", "object")
  # Replace numerical values with character labels
  long_data$Category <- category_labels[long_data$Category]
  
  # Convert input data from long format to wide format 
  wide_data <- pivot_wider(long_data, names_from = c(Species, Category), values_from = Measurement)
  
  # Save output 
  # Check if the file already exists
  if (file.exists("output.xlsx")) {
    # If the file exists, delete it
    file.remove("output.xlsx")
  }
  
  write.xlsx(wide_data, "output.xlsx", sheetName = "Input_data")
  
  # Summary statistics 
  summary_stats <- long_data %>%
    group_by(Species, Category) %>%
    get_summary_stats(Measurement, type = "mean_sd")
  write.xlsx(summary_stats, "output.xlsx", sheetName = "Summary_Stats", append = TRUE)
  
  
  # Visualize box plots 
  bxp <- ggboxplot(
    long_data, x = "Category", y = "Measurement",
    color = "Species", palette = "jco"
  )
  bxp
  
  # Check assumption #1: no outliers 
  # 'identify_outliers' function calculates Z-scores (# of SDs from the mean) for each data point, 
  # and data points > 3 SDs (or < -3 SDs) from the mean are considered as outliers 
  outliers <- long_data %>%
    group_by(Species, Category) %>%
    identify_outliers(Measurement)
  write.xlsx(outliers, "output.xlsx", sheetName = "Outliers", append = TRUE)
  
  # Check assumption #2: normality 
  # Shapiro wilk test: normal = p > 0.05
  normality <- long_data %>%
    group_by(Species, Category) %>%
    shapiro_test(Measurement)
  write.xlsx(normality, "output.xlsx", sheetName = "Normality", append = TRUE)
  
  
  # Perform rANOVA 
  results <- anova_test(
    data = long_data, dv = Measurement, wid = Subject,
    within = c(Species, Category)
  )
  
  # Show results 
  rmANOVA_results <- get_anova_table(results)
  write.xlsx(rmANOVA_results, "output.xlsx", sheetName = "rmANOVA_results", append = TRUE)
  
  
  # Post-hoc tests: Given IA, test for effect of species on category  
  pwc_results <- long_data %>%
    group_by(Category) %>%
    pairwise_t_test(
      Measurement ~ Species, paired = TRUE,
      p.adjust.method = "fdr"
    )
  pwc_results
  write.xlsx(pwc_results, "output.xlsx", sheetName = "pwc_results", append = TRUE)
  
  
  
  # Report results 
  pwc_results <- pwc_results %>% add_xy_position(x = "Category")
  
  bxp <- bxp + 
    stat_pvalue_manual(pwc_results, tip.length = 0, hide.ns = TRUE) +
    labs(
      caption = get_pwc_label(pwc_results)
    )
  
  # Save 
  # Check if the file already exists
  if (file.exists("boxplots.png")) {
    # If the file exists, delete it
    file.remove("boxplots.png")
  }
  ggsave("boxplots.png", plot = bxp, dpi = 300)
  
  

}
  

# Tutorial: https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/#two-way-repeated-measures-anova

