
rm(list=ls(all=TRUE))

# install any missing packages from CRAN
required_packages <- c("Matrix","TMB","fmesher","dplyr","mgcv","stats","tidyverse","data.table", 
                       "stringr","ggpubr","terra","tidyterra","RColorBrewer","patchwork","visreg")

missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}




#install specific sdmTMB version used (sensitive to version)

# unload the packages from memory for clean install
try(detach("package:sdmTMB", unload = TRUE), silent = TRUE)
try(unloadNamespace("sdmTMB"), silent = TRUE)

# install
unzipped_dir <- "./packages"

if (dir.exists(unzipped_dir)) {
  target_library <- .libPaths()[1]
  package_folders <- list.dirs(unzipped_dir, full.names = TRUE, recursive = FALSE)
  
  # get nested packages
  if (length(package_folders) == 1 && !(basename(package_folders) %in% c("sdmTMB"))) {
    package_folders <- list.dirs(package_folders, full.names = TRUE, recursive = FALSE)
  }
  
  # copies the package into your active R library
  for (folder in package_folders) {
    file.copy(from = folder, to = target_library, recursive = TRUE, overwrite = TRUE)
  }
}
