# download and install the ARTool package
# install.packages("ARTool")

# required packages 
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

# ART procedure is run on long-format data 
# load the long-format data and assign it to variable 'df'

setwd("/Users/jane_chesley/Documents/Github/Dataset_dynamic/EEG-analysis/statistics/R/PLI_normalized_roiHumanBody_Between-PLI_alpha")
file = "PLI_normalized_roiHumanBody_Between-PLI_alpha.xlsx"
df <- read_excel(file)

# Convert "Species" column to factors (strings)
df$Species <- factor(df$Species, levels = c(1, 2), labels = c("human", "monkey"))

# Convert "Category" column to factors (strings)
df$Category <- factor(df$Category, levels = c(1, 2, 3), labels = c("body", "face", "object"))


# 'Subject' is the name of a within-subjects column
# 'Species' is the name of the first factor column
# 'Category' is the name of the second factor column
# 'Measurement' is the name of the response column
# run the ART procedure on 'df'
m = art(Measurement ~ Species * Category + (1|Subject), data=df) # linear mixed-model syntax; see lme4::lmer()
anova(m)