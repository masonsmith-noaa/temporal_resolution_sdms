




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
#retrospective forecasting 
########################################################

rm(list=ls(pattern = "spearman_result.tmb."))
rm(list=ls(pattern = "rmse_result.tmb."))
rm(list=ls(pattern = "per_dev_expl_result.tmb."))

# load test data
region.data.test <- subset(region.data, year>2008)  

# run model loop
for(spp in top_spp) {
  for(step in c("annual", "timestep", "ltm", "ltm_static")){
    
    N.predictions_all = NULL
    N.model = readRDS(paste0("species_models/tmb_ar1_60_tw_m1/tmb.fit.ebs.", step, "_", spp, ".rds"))
    region.data.test$offset <- mean(N.model$offset)
    N.pred = predict(N.model, newdata = region.data.test, type="response") #counts
    N.obs.pred = N.pred
    N.obs.pred$forecast.yr = N.obs.pred$year
    N.obs.pred$fitted.through = 2023
    N.predictions_all = rbind(N.predictions_all, N.obs.pred)
    N.predictions_all$obs = N.predictions_all[,spp]
    N.predictions_all <- N.predictions_all %>% dplyr::select(x, y, year, est, obs) 
    
    # create region.data.preds for each year
    for(i in as.numeric(2009:2023)) {
      assign(paste0("region.data.preds_", i), N.predictions_all %>% dplyr::filter(year==i))
    }
    
    # assign 2020 for below loop to run, but don't include in final summary (no survey that year)
    region.data.preds_2020 <- region.data.preds_2019
    
    # Spearman's correlation
    for(i in 2009:2023) {
      assign(paste0("spearman_cor_", i), cor.test(formula = ~ obs + est, data = get(paste0("region.data.preds_", i)), method = c("spearman"), exact=FALSE))
    }
    
    
    spearman_cor_bind <- as.data.frame(rbind(spearman_cor_2009$estimate,spearman_cor_2010$estimate,
                                             spearman_cor_2011$estimate, spearman_cor_2012$estimate,
                                             spearman_cor_2013$estimate, spearman_cor_2014$estimate,
                                             spearman_cor_2015$estimate,spearman_cor_2016$estimate,
                                             spearman_cor_2017$estimate, spearman_cor_2018$estimate,
                                             spearman_cor_2019$estimate, spearman_cor_2021$estimate,
                                             spearman_cor_2022$estimate, spearman_cor_2023$estimate))
    
    spearman_cor_bind$year <- c(2009:2019, 2021:2023)
    spearman_cor_bind$spp <- spp
    spearman_cor_bind$step <- step
    
    assign(paste0("spearman_result.tmb.", step, ".", spp), spearman_cor_bind)
    
  }}   


# grab all results and bind into one dataframe and save
spearman_list <- mget(ls(pattern = "spearman_result.tmb."))
spearman_totals2.tmb_ar1_60_tw_m1 <- bind_rows(spearman_list)
saveRDS(spearman_totals2.tmb_ar1_60_tw_m1, "data/spearman_totals2.tmb_ar1_60_tw_m1.rds")