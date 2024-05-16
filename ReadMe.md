<!--- This is how you write comments that do not appear to reader --->

<!-- back to top link -->
<a name="readme-top"></a>

# EEG-analysis 



<!-- TABLE OF CONTENTS -->
### Table of Contents 
  <ol>
    <li><a href="#description">Description</a></li>
    <li><a href="#getting-started">Getting Started</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>

<br>

<!-- ABOUT THE PROJECT -->
## Description  

This repository provides a comprehensive pipeline for the analysis of 5-D electroencephalography (EEG) data. The pipeline includes preprocessing, spectral analyis, connectivity analysis, statistical analysis and data visualization. Analysis of EEG data can reveal meaningful insights of brain activity.  

<p align="right">(<a href="#readme-top">back to top</a>)</p>
<br>

<!-- Getting started  -->
## Getting Started

To get started with this pipeline, follow these steps:
1. Install MATLAB on your system, along with the following toolboxes:
    - Image Processing Toolbox
    - Statistics and Machine Learning Toolbox
    - Bioinformatics Toolbox
    - FieldTrip Toolbox
2. Install RStudio on your system, along with the following packages:
    - ggplot2
    - ggpubr
    - rstatix
    - dplyr
    - tidyverse
    - xlsx
    - readxl
    - knitr
    - broom
    - car 
3. Download or clone this repository, 'EEG-analysis'
4. Launch MATLAB and navigate to the repository directory '/main_scripts'
5. Save raw EEG data in the repository directory '/input' 
6. Preprocess raw EEG data with the script 'preprocessing.m', which includes:
    - Trial segmentation and extraction
    - Event labeling
    - Signal re-referencing
    - Sampling reduction
    - Electrooculogram (EOG) artifact removal 
    - Removal of outlier trials
    - Stratification by condition 
7. Perform time-frequency power analysis of preprocessed EEG data with the script 'TF_power_analysis.m', which includes:
    - Bandpass filtering
    - Morlet wavelet transformation 
    - Baseline correction  
    - Data visualizations
8. Perform phase-based connectivity analysis with the script 'connectivity_analysis.m'. which includes:
    - Single-trial 5-D Phase-Lag Index (PLI) computations 
    - PLI computation within-regions
    - PLI computation between-regions
    - Validation of computations with simulated data 
9. View all processed data in the repository directory '/output'
8. Launch RStudio and navigate to the repository directory '/statistics' 
9. Execute 'Stats.R', which includes:
    - Repeated measures ANOVA on PLI data 
    - Post-hoc paired t-tests 
    - Data visualizations 
10. View all statistical results and plots in the repository directory '/statistics/R' 

<p align="right">(<a href="#readme-top">back to top</a>)</p>
<br>

<!-- VISUALIZE PROJECT OUTLINE  -->
## Roadmap

Visualizations of project outline.

<p align="right">(<a href="#readme-top">back to top</a>)</p>
<br>

<!-- USAGE  -->
## Usage 

This usage of this pipeline is two-fold. First, the pipeline is a comprehensive toolkit for understanding and analyzing multidimensional EEG data. Second, the pipeline can be used to replicate results of the study: 'Dynamic Body Processes in the Human Brain.' To replicate results of that study, please reach out to the author to initiate a formal data sharing agreement. 

<p align="right">(<a href="#readme-top">back to top</a>)</p>
<br>

<!-- CONTRIBUTING -->
## Contributing

Contributions are appreciated. If you have a suggestion, please fork the repo and create a pull request or reach out to the author. 

<p align="right">(<a href="#readme-top">back to top</a>)</p>
<br>


<!-- CONTACT -->
## Contact

### Your Name <br>
[janechesley@gmail.com](janechesley@gmail.com) <br>
[LinkedIn](https://www.linkedin.com/in/jane-chesley/) <br> <br>
### Project Link: 
[https://github.com/Jane-Chesley/EEG-analysis.git](https://github.com/Jane-Chesley/EEG-analysis.git)

<p align="right">(<a href="#readme-top">back to top</a>)</p>
<br>