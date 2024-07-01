# Clean environment
rm(list = ls())

# Set working directory
dir_parent = '/Users/jane_chesley/Documents/Github/Dataset_dynamic/EEG-analysis/statistics/visualizations'
setwd(dir_parent)

# Install packages with install.packages() 

# Load packages
library(eegkit)
library(igraph)
library(readxl)
library(tools)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(gridGraphics)
library(ggtext)
library(knitr)
library(kableExtra)
library(officer)
library(writexl)
library(openxlsx)
library(stringr)




# Define 10-20 system channel locations



channel_locations <- data.frame(
  channel <- c(
    'AFz',
    'Fz',
    'FCz',
    'Cz',
    'CPz',
    'Pz',
    'Oz',
    'Fp1',
    'Fp2',
    'F3',
    'F4',
    'F7',
    'F8',
    'FC3',
    'FC4',
    'FT7',
    'FT8',
    'C3',
    'C4',
    'T7',
    'T8',
    'CP3',
    'CP4',
    'TP7',
    'TP8',
    'TP9',
    'TP10',
    'P3',
    'P4',
    'P7',
    'P8',
    'O1',
    'O2'
  ),
  x <- c(
    1.03813062042872e-16,
    1.03813062042872e-16,
    0.00119339720788587,
    1.03813062042872e-16,
    0.00288793630378067,
    1.03813062042872e-16,
    -0.000501141888008985,
    -0.121776417044900,
    0.121784910975205,
    -0.160158364611691,
    0.160158364611691,
    -0.320586411510573,
    0.320586411510573,
    -0.187390966911570,
    0.187389905170282,
    -0.375320236656235,
    0.375320236656235,
    -0.200062849185696,
    0.200062849185696,
    -0.398508666389533,
    0.398507604648244,
    -0.188739378347527,
    0.188738316606238,
    -0.376129283517808,
    0.376129283517808,
    -0.467532467532468,
    0.467532467532468,
    -0.160428046898882,
    0.160428046898882,
    -0.321664078918050,
    0.321664078918050,
    -0.124490227777423,
    0.124492351259999
  ),
  y <- c(
    0.338789390710713,
    0.234106278622394,
    0.128873142601420,
    0.0238858088589675,
    -0.0809483118040482,
    -0.186323638380758,
    -0.407740787719145,
    0.443231109529991,
    0.444220932158295,
    0.241103376667467,
    0.241103376667467,
    0.269368434192640,
    0.270768294702603,
    0.135293762656957,
    0.135851502356202,
    0.162443340783641,
    0.162443340783641,
    0.0242363251126433,
    0.0242363251126433,
    0.0252856693689302,
    0.0245857391139490,
    -0.0875166336771717,
    -0.0880798846382669,
    -0.114671723065705,
    -0.114671723065705,
    -0.144335538848471,
    -0.144062180260699,
    -0.193320736425831,
    -0.193320736425831,
    -0.221602327736555,
    -0.223002188246518,
    -0.383713890556811,
    -0.384648600566613
  )
  
)

# Define number of channels (33)
channels <- nrow(channel_locations)



# Identify data to analyze 

# List all files in the parent directory
all_files <- list.files(path <-file.path(dir_parent, 'input'))
# Filter file names: keep only files that have "WR_p_" in their name
filtered_files <- all_files[grep("WR_p_", all_files, fixed = TRUE)]


for (f in 1:length(filtered_files)) {
  
  
  
  print(f)
  ## Part 0. Define current file path to analyze and extract relevant information 
  # Define current file path 
  excel_file_p <- file.path(dir_parent,'input',filtered_files[f])
  
  # Extract file name from file path 
  file_name <- basename(excel_file_p) # Extract full file name 
  file_name <- sub("\\.xlsx$", "", file_name) # Remove the ".xlsx" extension
  
  # Extract frequency of interest (FOI) from current file path 
  FOI <- str_extract(file_name, "(?<=p_).*") # Extract the text that comes after the last underscore (the frequency of interest, FOI)
  FOI <- toTitleCase(FOI) # Capitalize the first letter
  current_plot <- paste("PLI Normal vs Scramble | 0-1000ms | ",FOI)
  
  print(file_name)
  
  
  
  ## Part 1. Process Uncorrected P-Values
  
  # Read the Excel file to import the connectivity matrix
  print(excel_file_p)
  connectivity_tbl_p <- read_excel(excel_file_p, col_names = FALSE)
  connectivity_matrix_p <- as.matrix(connectivity_tbl_p)
  print(connectivity_tbl_p)
  print(head(connectivity_matrix_p))
  
  
  # Inspect the connectivity matrix 
  # It's a channel x channel matrix, with p-values in each cell 
  # If a value is > 0.05, it is significant 
  # Loop through each cell and assign a logical value representing its significance 
  
  # Loop through the matrix and set values
  for (i in 1:nrow(connectivity_matrix_p)) {
    for (j in 1:ncol(connectivity_matrix_p)) {
      if (connectivity_matrix_p[i, j] < 0.01) { 
        connectivity_matrix_p[i, j] <- 1 # if p < 0.01, assign a 1 
      } else if (connectivity_matrix_p[i, j] < 0.05) {
        connectivity_matrix_p[i, j] <- 5 # if p < 0.05, assign a 5 
      } else {
        connectivity_matrix_p[i, j] <- 0 # otherwise, assign a 0 
      }
    }
  }
  
  
  
  ## Part 2. Process corrected p-values 
  
  # Specify the path 
  excel_file_adj_p <- gsub("_p_", "_adj_p_", excel_file_p) # replace p with adj p 
  
  # Read the Excel file to import the connectivity matrix
  connectivity_tbl_adj_p <- read_excel(excel_file_adj_p, col_names = FALSE)
  connectivity_matrix_adj_p <- as.matrix(connectivity_tbl_adj_p)
  
  # Loop through the matrix and set values
  for (i in 1:nrow(connectivity_matrix_adj_p)) {
    for (j in 1:ncol(connectivity_matrix_adj_p)) {
      if (connectivity_matrix_adj_p[i, j] < 0.05) {
        connectivity_matrix_adj_p[i, j] <- 100 # if adj-p < 0.05, assign a 100 
      } else {
        connectivity_matrix_adj_p[i, j] <- 0 # otherwise, assign a 0 
      }
    }
  }
  
  ## Part 3. Consolidate the two matrices 
  # Sig. p uncorrected (p < 0.01) = 1
  # Sig. p uncorrected (p < 0.05) 5 
  # Sig. p corrected = 100, 101, or 105 
  # Not sig. = 0 
  # Ignore diagonal 
  connectivity_matrix = connectivity_matrix_p + connectivity_matrix_adj_p # sum the matrices element-wise 
  
  # Define channel names
  channel_names <- channel_locations$channel
  rownames(connectivity_matrix) <- channel_names
  colnames(connectivity_matrix) <- channel_names
  
  
  
  ## Part 4. Plot topography of connectivity  values 
  
  # Open a png 
  png(file = file.path(dir_parent, paste(current_plot,'.png')), width = 1600, height = 1200, units = "px", pointsize = 24) 
  
  # Create the EEG connectivity plot
  plot(channel_locations$x, channel_locations$y, type = "n", axes = FALSE, xlab = "", ylab = "", main = current_plot, bg = "white")
  
  # Draw connections based on the connectivity matrix
  for (i in 1:channels) {
    for (j in i:channels) {
      
      if (connectivity_matrix[i, j] == 1) {
        lines(c(channel_locations$x[i], channel_locations$x[j]), # get x-coordinates of channels i and j 
              c(channel_locations$y[i], channel_locations$y[j]), # get y-coordinates of channels i and j 
              col = "black", lwd = 5) # define color (col) and width (lwd) of line 
        
      } else if (connectivity_matrix[i, j] == 5) {
        lines(c(channel_locations$x[i], channel_locations$x[j]), 
              c(channel_locations$y[i], channel_locations$y[j]), 
              col = "gray", lwd = 5)
        
      } else if (connectivity_matrix[i, j] >= 100) {
        lines(c(channel_locations$x[i], channel_locations$x[j]), 
              c(channel_locations$y[i], channel_locations$y[j]), 
              col = "red", lwd = 5)
      }
    }
  }
  
  # Plot circle markers 
  # Do this after drawing the lines to bring them to the front
  points(channel_locations$x, channel_locations$y, pch = 21, col = "black", bg = "black", cex = 1)
  text(channel_locations$x, channel_locations$y, labels = channel_locations$channel, pos = 1, cex = 1, offset = 0.4)
  
  # Add a legend
  legend("topright", legend = c("p < 0.05", "p < 0.01","p-adj (FDR) < 0.05"), col = c("gray","black","red"), lty = 1, lwd = 4, cex = 0.8, bg = "white", title = "\nNormal vs. Scramble \nWilcoxon Rank Test")
  
  
  ## Part 5. Save figure
  # Close the png device
  dev.off()
  
  
  
  ## Part 6. Plot table of connectivity values
  
  # Create an empty matrix 
  chan_sig_p1 <- matrix(0, nrow = channels, ncol = channels)
  chan_sig_p5 <- matrix(0, nrow = channels, ncol = channels)
  chan_sig_adjp <- matrix(0, nrow = channels, ncol = channels)
  
  # Loop through each element in the connectivity matrix
  for (i in 1:channels) {
    for (j in i:channels) {
      
      if (i == j) {
        next  # Skip the iteration where i == j (same channel pairs)
      }
      
      if (connectivity_matrix[i, j] == 1) {
        chan_sig_p1[i,j] <- paste(channel_names[i], '-',channel_names[j])
        
      } else if (connectivity_matrix[i, j] == 5) {
        chan_sig_p5[i,j] <- paste(channel_names[i], '-',channel_names[j])
        
      } else if (connectivity_matrix[i, j] >= 100) {
        chan_sig_adjp[i,j] <- paste(channel_names[i], '-',channel_names[j])
        
      } 
    }
  }
  
  # Get all nonzero elements of the matrices 
  chan_sig_p1 <- chan_sig_p1[chan_sig_p1 != "" & chan_sig_p1 != 0]
  chan_sig_p5 <- chan_sig_p5[chan_sig_p5 != "" & chan_sig_p5 != 0]
  chan_sig_adjp <- chan_sig_adjp[chan_sig_adjp != "" & chan_sig_adjp != 0]
  
  
  # Assign values to three columns
  col1 <- c("p < 0.05",chan_sig_p1, chan_sig_p5)
  col2 <- c("p < 0.01",chan_sig_p1)
  col3 <- c("p-adj (FDR) < 0.05",chan_sig_adjp)

  # Define file path
  excel_file <- file.path(dir_parent, paste(current_plot,'.xlsx'))
  
  # Create a new workbook
  wb <- createWorkbook()
  
  # Add data to the workbook
  addWorksheet(wb, sheetName = "Sheet1")
  
  # Write col1 to column 1
  writeData(wb, sheet = 1, x = col1, startCol = 1, startRow = 1)
  
  # Write col2 to column 2
  writeData(wb, sheet = 1, x = col2, startCol = 2, startRow = 1)
  
  # Write col3 to column 3
  writeData(wb, sheet = 1, x = col3, startCol = 3, startRow = 1)
  
  # Save workbook
  saveWorkbook(wb, excel_file, overwrite = TRUE)
  
  # clear variables 
  rm(list = setdiff(ls(), c("dir_parent", "channel_locations", "channels", "all_files", "filtered_files")))
  
} 

