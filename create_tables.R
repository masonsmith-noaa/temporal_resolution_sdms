rm(list=ls(all=TRUE))

# install any missing required packages and load libraries
required_packages <- c("sdmTMB", "mgcv", "tidyverse", "data.table")
lapply(required_packages, library, character.only = TRUE)

#set parameters
region.name = "EBS"

#load catch data ==============================================================================
region.data <- read.csv("data/region.data.csv", sep=",", header=TRUE)
region.data$sponge <- as.factor(region.data$sponge)
region.data$pen <- as.factor(region.data$pen)
region.data$coral <- as.factor(region.data$coral)

########################################################
#filter top frequently occurring species 
########################################################

# count occurrences
region.data.freq <- as.data.frame(colSums(region.data[, 52:161] > 0))
region.data.freq2 <- cbind(rownames(region.data.freq), data.frame(region.data.freq, row.names = NULL))
colnames(region.data.freq2) <- c("spp", "frequency")
region.data.freq2 <- subset(region.data.freq2, spp != "halibut") 

# arrange by frequency and select the top 20
region.data.freq3 <- region.data.freq2 %>%
  dplyr::arrange(desc(frequency)) %>%
  dplyr::slice(1:15)

top_spp <- paste0(region.data.freq3$spp)
print(top_spp)


#add full names
spp_names_df <- data.frame(
  species = top_spp,
  full_name = c("Adult Pollock","Juv. Pacific Cod","Adult Pacific Cod","Snow Crab","Adult Flathead Sole",
                "Juv. Flathead Sole","Juv. Pollock","Early Juv. Pollock","Tanner Crab","Adult Yellofin Sole",
                "Juv. Yellowfin Sole","Adult Alaska Plaice","Adult N. Rock Sole","Juv. N. Rock Sole","Juv. Alaska Skate")
)



#########################################################################
### Table S1 - summarizes range estimates from Matern functions 
#########################################################################

# empty dataframe to store results
all_ran_pars <- data.frame()

# loop through species
for(spp in top_spp) {
  tryCatch({
    file_path <- paste0("species_models/tmb_ar1_60_tw_m1/tmb.fit.ebs.null_", spp, ".rds")
    N.model <- readRDS(file_path)
    ran_pars <- tidy(N.model, "ran_pars") #random effects
    ran_pars$species <- spp
    ran_pars$model <- "null"
    
    #combine to main dataframe 
    all_ran_pars <- rbind(all_ran_pars, ran_pars)
    
  }, error = function(e) {
    message(paste("Skipping model for species:", spp, "due to error:", e$message))
  })
}

#save ranges
ranges <- all_ran_pars %>% filter(term == "range")
ranges <- ranges %>% dplyr::left_join(spp_names_df, by="species")
write.csv(ranges, paste0("tables/Table_S1_ranges.csv"))



###############################################################
### Table S2 - Sanity Check Results
###############################################################

#load models and run sanity() checks  
sanity_results <- list() 

  for(spp in top_spp) {
    
    sanity_results[[spp]] <- list()
    
    #annual----------------------------------------------------------------------------
    tmb.fit.ebs.annual <- readRDS(paste0("species_models/tmb_ar1_60_tw_m1/tmb.fit.ebs.annual_", spp, ".rds"))
    print("Running sanity check for annual model...")
    sanity_results[[spp]]$annual <- sanity(tmb.fit.ebs.annual)

    #5-year timestep-------------------------------------------------------------------
    tmb.fit.ebs.timestep <- readRDS(paste0("species_models/tmb_ar1_60_tw_m1/tmb.fit.ebs.timestep_", spp, ".rds"))
    print("Running sanity check for 5-year timestep model...")
    sanity_results[[spp]]$timestep <- sanity(tmb.fit.ebs.timestep)

    #LTM-------------------------------------------------------------------------------
    tmb.fit.ebs.ltm <- readRDS(paste0("species_models/tmb_ar1_60_tw_m1/tmb.fit.ebs.ltm_", spp, ".rds"))
    print("Running sanity check for LTM model...")
    sanity_results[[spp]]$ltm <- sanity(tmb.fit.ebs.ltm)
    
    #LTM_static-------------------------------------------------------------------------------
    tmb.fit.ltm_static <- readRDS(paste0("species_models/tmb_ar1_60_tw_m1/tmb.fit.ebs.ltm_static_", spp, ".rds"))
    print("Running sanity check for LTM_static model...")
    sanity_results[[spp]]$ltm_static <- sanity(tmb.fit.ltm_static)
    
  } 
  
  
#empty dataframe to store results
sanity_df <- data.frame(
  species = character(),
  model = character(),
  check = character(),
  result = logical()
)

# loop through species results
for (spp in names(sanity_results)) {
  for (model_name in names(sanity_results[[spp]])) {
    current_results <- sanity_results[[spp]][[model_name]]
    if (!is.null(current_results)) {
      checks <- names(current_results)
      results <- unlist(current_results)
      
      temp_df <- data.frame(
        species = spp,
        model = model_name,
        check = checks,
        result = results
      )
      
      sanity_df <- rbind(sanity_df, temp_df)
    }
  }
}

# save results
sanity_df_wide <- sanity_df %>%
  pivot_wider(
    names_from = "check",
    values_from = "result"
  )

print("Sanity Check:")
print(sanity_df_wide)

write.csv(sanity_df_wide, "tables/Table_S1_Sanity_Check.csv")
  


###############################################################
### Table S3 - Model Fit Summaries
###############################################################

# define models and target random parameters
model_steps <- c("annual", "timestep", "ltm", "ltm_static")
target_terms <- c("sigma_O", "sigma_E", "sigma_Z", "ar1_rho")

# loop
all_extracted_data <- data.frame()

for (step in model_steps) {
  for (spp in top_spp) {
    file_path <- paste0("species_models/tmb_ar1_60_tw_m1/tmb.fit.ebs.", step, "_", spp, ".rds")
    if (!file.exists(file_path)) next
    fit <- readRDS(file_path)
    
    # extract fixed effects and smooths
    model_report <- fit$sd_report$par.fixed
    report_df <- data.frame(term = names(model_report), estimate = model_report, stringsAsFactors = FALSE)
    smooth_df <- report_df[grep("^ln_smooth_sigma", report_df$term), ] 
    
    if (nrow(smooth_df) == 3) {
      smooth_df$estimate <- exp(smooth_df$estimate) 
      smooth_df$term <- c("sds(bdepth_scaled)", "sds(bcurrentU_scaled,bcurrentV_scaled)", "sds(btemp_scaled)")
      smooth_df$measure <- "SD"
      smooth_df$description <- "effect of smooth"
    } else {
      smooth_df <- data.frame(term = paste0("smooth_", 1:nrow(smooth_df)), estimate = exp(smooth_df$estimate), measure = "SD", description = "effect of smooth")
    }
    
    # random effects
    st_raw <- as.data.frame(tidy(fit, "ran_pars"))
    st_df <- st_raw[st_raw$term %in% target_terms, c("term", "estimate")]
    
    if (nrow(st_df) > 0) {
      st_df$measure <- ifelse(st_df$term == "ar1_rho", "rho", "SD")
      st_df$description <- case_when(
        st_df$term == "sigma_O" ~ "spatial random effect",
        st_df$term == "sigma_E" ~ "AR1 random effect",
        st_df$term == "sigma_Z" ~ "variation of SVC",
        st_df$term == "ar1_rho" ~ "temporal correlation"
      )
    }
    
    # combine all models
    combined_df <- rbind(smooth_df, st_df)
    combined_df$step <- step
    combined_df$species <- spp
    all_extracted_data <- rbind(all_extracted_data, combined_df)
  }
}

# join names
final_flat_df <- all_extracted_data %>% dplyr::left_join(spp_names_df, by = "species")

# pivot wider by timestep
final_total_df_wide <- final_flat_df %>%
  pivot_wider(
    names_from = step,
    values_from = estimate)

# save final table (Table S3)
print(final_total_df_wide)
write.csv(final_total_df_wide, paste0("tables/Table_S3_Model_Fits.csv"))
