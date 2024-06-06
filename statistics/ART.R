# this script runs ART on PLI datasets 

# ------------------------------------------------
# STEP 1: SETUP ENVIRONMENT
# ------------------------------------------------

# clean environment 
rm(list = ls())

# load required packages 
library(ARTool)
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

# define parent directory
dir_parent = "/Users/jane_chesley/Documents/Github/Dataset_dynamic/EEG-analysis/statistics/R"



# ------------------------------------------------
# STEP 2: PREPARE DATA 
# ------------------------------------------------

# list all folders in the parent directory
# do not list subfolders ('recursive = FALSE')
all_folders <- list.dirs(dir_parent, full.names = FALSE, recursive = FALSE)

# filter folder names: keep only folders that have "PLI_normalized" in their name
filtered_folders <- all_folders[grep("PLI_normalized", all_folders, fixed = TRUE)]

# construct file names from folder names 
file_names <- paste(filtered_folders, ".xlsx", sep = "")

# identify file names and folders of all data to analyze 
all_data <- data.frame(folder_name = filtered_folders, file_name = file_names)


# ------------------------------------------------
# STEP 3: ART PROCEDURE  
# ------------------------------------------------

# start loop to analyze all data  
i <-  1
current_folder <- all_data[i, "folder_name"]
current_file <- all_data[i, "file_name"]

# set working directory  
setwd(paste(dir_parent, current_folder, sep = "/"))

# ART procedure requires:
    # long-format data  
    # factors (strings) in subject and condition columns 
    # first column = subject identifier 'S'
    # last column = numeric response variable 'Y'
    # all columns between first and last are independent variables (IV), E.g. X1, X2, X3

# load input data in long format and assign it to df 
df <- read_excel(current_file)

# inspect the format of df  
str(df)
head(df, n=10) # display only the first 10 observations

# the present datasets have numbers (format <dbl>) in subject and condition columns, which is not suitable for ART 
# convert the columns to factors (strings; format <fct>)
df$Subject <- factor(df$Subject)
df$Species <- factor(df$Species, levels = c(1, 2), labels = c("human", "monkey"))
df$Category <- factor(df$Category, levels = c(1, 2, 3), labels = c("body", "face", "object"))

# inspect the new format of df  
str(df)
head(df, n=10) # display only the first 10 observations

# variable 'df' contains: 
    # Col1: Subject 
    # Col2: Species (IV 1), with levels 'human' and 'monkey
    # Col3: Category (IV 2), with levels 'body', 'face', and 'object'
    # Col4: 'Measurement', which represents the response variable of PLI values

# transform the data, to be analyzed with a repeated measures ANOVA (rmANOVA)
# rmANOVA includes an Error term 
transformedDataRM <- art(Measurement ~ Species*Category + Error(Subject), data=df)

# verify the transformation was done correctly:
    # the aligned responses should = 0 
    # the F values of ANOVAS on aligned responses should = 0 
    # zero values mean the alignment correctly "stripped out" effects not of interest 

# inspect the transformed data 
summary(transformedDataRM)

# if the data was transformed correctly, proceed with the statistical test, rmANOVA 
anova(transformedDataRM)

# transform the data, to be analyzed with a mixed effects model 
transformedDataME <- art(Measurement ~ Species*Category + (1|Subject), data=df)

# inspect the transformed data 
summary(transformedDataME)

# if the data was transformed correctly, proceed with the statistical test, mixed effects model 
anova(transformedDataME)
