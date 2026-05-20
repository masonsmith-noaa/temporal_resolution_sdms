This repository includes the data and code for running the scripts in: "Effects of temporal resolution on forecasting species distributions under varying conditions." To run all scripts, unzip this folder to a permanent location and create a new R Project out of this folder. In RStudio, this can be done through Project > New Project > Existing Directory > Browse to Folder > Open. The code will now draw from all folders within the source folder.

Files and variables

File: data folder

Description: This folder contains the source file containing catch and environmental data ('region.data.csv') to run each spatiotemporal model, and the predictive grid ('grid/grid.agg10.newbathy.ebs.rds') for producing raster maps. This folder also serves as storage for hosting predictions ('predictions' subfolder) and retrospective forecast results ('spearman_totals2.tmb_ar1_60_tw_m1.rds').

File: figures folder

Description: This folder serves as storage for all figures saved from 'create_figures.R'.

File: tables folder

Description: This folder serves as storage for all tables saved from 'create_tables.R' and the sdmTMB sanity check from 'model_runs.R'.

File: species_models folder

Description: This folder serves as storage for all models saved in 'model_runs.R'.

File: model_runs.R

Description: This script runs the spatiotemporal models for each temporal resolution (annual, timestep, long-term mean, and static). Each model will be saved under the 'species_models' subfolder. Sanity checks from sdmTMB will additionally be saved under the 'sanity_tables' subfolder.

File: retrospective_forecasts.R

Description: This script runs the retrospective forecasts for each temporal resolution (annual, timestep, long-term mean, and static) and saves a final dataframe under the 'data' subfolder.

File: create_tables.R

Description: This script creates tables S1 (Matérn range estimates), S2 (sanity check results), and S3 (model fit summaries) of the supplemental appendix.

File: create_figures.R

Description: This script creates figures 3 (maximum Cohen's q comparisons), 4 (Spearman's rho of retrospective forecasts), 5 (Pairwise differences in Fisher’s-standardized forecast skill), 6 (Spatial agreement of core habitat predictions among temporal resolutions), and 7 (Predicted spatial agreement of core habitat) of the main text; and figures S1 (spatial mesh configuration), S2 (Conditional response plots), S3 (Spearman’s rho of retrospective forecasts), and S4 (Predicted area of 50th percentile in abundance among the four temporal resolutions) of the supplemental appendix.

Code/software

R 4.4.3

RStudio 2024.12.1

required R packages: "Matrix", "TMB", "sdmTMB", "fmesher", "sdmTMBextra", "dplyr", "mgcv", "stats", "tidyverse", "data.table", "stringr", "sdmTMB", "ggpubr", "terra", "tidyterra","RColorBrewer", "tidyverse", "data.table", "patchwork", "visreg"



