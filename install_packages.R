
rm(list=ls(all=TRUE))

# install any missing required packages and load libraries
required_packages <- c("sdmTMB","Matrix","TMB","fmesher","dplyr","mgcv","stats","tidyverse","data.table", 
                       "stringr","ggpubr","terra","tidyterra","RColorBrewer","patchwork","visreg")

missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}