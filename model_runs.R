

rm(list=ls(all=TRUE))

# install any missing required packages and load libraries
required_packages <- c(
  "Matrix", "TMB", "sdmTMB", "fmesher", "sdmTMBextra", 
  "dplyr", "mgcv", "stats", "tidyverse", "data.table", "stringr"
)

missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

lapply(required_packages, library, character.only = TRUE)


#set parameters
region.name = "EBS"


#load catch data ==============================================================================
region.data <- read.csv("data/region.data.csv", sep=",", header=TRUE)


########################################################
#filter top frequently occurring species 
########################################################

# 1. Count the number of non-zero occurrences for each species
region.data.freq <- as.data.frame(colSums(region.data[, 52:161] > 0))
region.data.freq2 <- cbind(rownames(region.data.freq), data.frame(region.data.freq, row.names = NULL))
colnames(region.data.freq2) <- c("spp", "frequency")
region.data.freq2 <- subset(region.data.freq2, spp != "halibut")

# Arrange by frequency and select the top 20
region.data.freq3 <- region.data.freq2 %>%
  dplyr::arrange(desc(frequency)) %>%
  dplyr::slice(1:15) 

top_spp <- paste0(region.data.freq3$spp)
print(top_spp)

train.data = subset(region.data, year < 2009)

########################################################
#make mesh for tmb models and cross-validation
########################################################

mesh <- sdmTMB::make_mesh(train.data, c("x", "y"),
                          fmesher_func = fmesher::fm_mesh_2d_inla,
                          cutoff = 0.6, # minimum triangle edge length
                          max.edge = c(3, 9),
                          offset=c(1.2, 9)) # inner and outer max triangle lengths
plot(mesh)
mesh$mesh$n 



########################################################
#run models 
########################################################
for(spp in top_spp) {
  
  sanity_results[[spp]] <- list() 
  
  #null----------------------------------------------------------------------------
  tmb.form.null <- formula(paste0(spp, " ~ 1"))

  tmb.fit.ebs.null <- sdmTMB::sdmTMB(
    data = train.data,
    formula = tmb.form.null,
    family = tweedie(link="log"),
    spatial = "off",
    offset = "logarea"
  )

  saveRDS(tmb.fit.ebs.null, paste0("species_models/tmb_ar1_60_tw_m1/tmb.fit.ebs.null_", spp, ".rds"))


  
  #annual----------------------------------------------------------------------------

  tmb.form.annual <- formula(paste0(spp, " ~ s(bdepth_scaled, bs='tp', m=1, k=7) +
                                               s(bcurrentU_scaled, bcurrentV_scaled, bs='tp', m=1, k=7) +
                                               s(btemp_scaled, bs='tp', m=1, k=7) +
                                               cpe_scaled"))

  tmb.fit.ebs.annual <- sdmTMB::sdmTMB(
    data = train.data,
    formula = tmb.form.annual,
    mesh = mesh,
    family = tweedie(link="log"),
    spatial = "on",
    time = "year",
    extra_time = 2009:2023,
    spatial_varying = ~ 0 + cpe_scaled,
    offset = "logarea",
    spatiotemporal = "ar1"
  )

  saveRDS(tmb.fit.ebs.annual, paste0("species_models/tmb_ar1_60_tw_m1/tmb.fit.ebs.annual_", spp, ".rds"))


  # 5-year timestep-------------------------------------------------------------------
  tmb.form.timestep <- formula(paste0(spp, " ~ s(bdepth_scaled, bs='tp', m=1, k=7) +
                                               s(bcurrentU_timestep_scaled, bcurrentV_timestep_scaled, bs='tp', m=1, k=7) +
                                               s(btemp_timestep_scaled, bs='tp', m=1, k=7) +
                                               cpe_timestep_scaled"))

  tmb.fit.ebs.timestep <- sdmTMB::sdmTMB(
    data = train.data,
    formula = tmb.form.timestep,
    mesh = mesh,
    family = tweedie(link="log"),
    spatial = "on",
    time = "year",
    extra_time = 2009:2023,
    spatial_varying = ~ 0 + cpe_timestep_scaled,
    offset = "logarea",
    spatiotemporal = "ar1"
  )

  saveRDS(tmb.fit.ebs.timestep, paste0("species_models/tmb_ar1_60_tw_m1/tmb.fit.ebs.timestep_", spp, ".rds"))



  # LTM-------------------------------------------------------------------------------
  tmb.form.ltm <- formula(paste0(spp, " ~  s(bdepth_scaled, bs='tp', m=1, k=7) +
                                             s(bcurrentU_ltm_scaled, bcurrentV_ltm_scaled, bs='tp', m=1, k=7) +
                                             s(btemp_ltm_scaled, bs='tp', m=1, k=7)"))

  tmb.fit.ebs.ltm <- sdmTMB::sdmTMB(
    data = train.data,
    formula = tmb.form.ltm,
    mesh = mesh,
    family = tweedie(link="log"),
    spatial = "on",
    time = "year",
    extra_time = 2009:2023,
    offset = "logarea",
    spatiotemporal = "ar1"
  )

  saveRDS(tmb.fit.ebs.ltm, paste0("species_models/tmb_ar1_60_tw_m1/tmb.fit.ebs.ltm_", spp, ".rds"))



  # LTM_static-------------------------------------------------------------------------------
  tmb.form.ltm_static <- formula(paste0(spp, " ~  s(bdepth_scaled, bs='tp', m=1, k=7) +
                                           s(bcurrentU_ltm_scaled, bcurrentV_ltm_scaled, bs='tp', m=1, k=7) +
                                           s(btemp_ltm_scaled, bs='tp', m=1, k=7)"))

  tmb.fit.ltm_static <- sdmTMB::sdmTMB(
    data = train.data,
    formula = tmb.form.ltm_static,
    mesh = mesh,
    family = tweedie(link="log"),
    spatial = "on",
    offset = "logarea"
  )

  saveRDS(tmb.fit.ltm_static, paste0("species_models/tmb_ar1_60_tw_m1/tmb.fit.ebs.ltm_static_", spp, ".rds"))
  
}







