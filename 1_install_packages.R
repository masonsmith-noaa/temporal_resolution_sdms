
rm(list=ls(all=TRUE))

# install any missing packages from CRAN
required_packages <- c("Matrix","TMB","fmesher","dplyr","mgcv","stats","tidyverse","data.table", 
                       "stringr","ggpubr","terra","tidyterra","RColorBrewer","patchwork","visreg")

missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}



#install specific sdmTMB version used (sensitive to version)

# 1. Unload the packages from memory so Windows releases the .dll files
try(detach("package:sdmTMB", unload = TRUE), silent = TRUE)
try(unloadNamespace("sdmTMB"), silent = TRUE)

# 2. Run the copy script again
zip_file <- "./packages.zip"
temp_unzip_dir <- "./temp_packages_extracted"

if (file.exists(zip_file)) {
  unzip(zipfile = zip_file, exdir = temp_unzip_dir)
  target_library <- .libPaths()[1]
  package_folders <- list.dirs(temp_unzip_dir, full.names = TRUE, recursive = FALSE)
  
  if (length(package_folders) == 1 && !(basename(package_folders) %in% c("sdmTMB"))) {
    package_folders <- list.dirs(package_folders, full.names = TRUE, recursive = FALSE)
  }
  
  for (folder in package_folders) {
    file.copy(from = folder, to = target_library, recursive = TRUE, overwrite = TRUE)
  }
  
  unlink(temp_unzip_dir, recursive = TRUE)
}
