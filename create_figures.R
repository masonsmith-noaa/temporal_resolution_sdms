rm(list=ls(all=TRUE))

#libraries
required_packages <- c("sdmTMB", "ggpubr", "mgcv", "stats", "terra", "tidyterra", 
  "RColorBrewer", "tidyverse", "data.table", "patchwork", "visreg")

missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)}
lapply(required_packages, library, character.only = TRUE)

#set region
region.name = "EBS"

#load catch data ==============================================================================
region.data <- read.csv("data/region.data.csv", sep=",", header=TRUE)
region.data$sponge <- as.factor(region.data$sponge)
region.data$pen <- as.factor(region.data$pen)
region.data$coral <- as.factor(region.data$coral)

########################################################
#filter top X species 
########################################################

# count occurrences
region.data.freq <- as.data.frame(colSums(region.data[, 52:161] > 0))
region.data.freq2 <- cbind(rownames(region.data.freq), data.frame(region.data.freq, row.names = NULL))
colnames(region.data.freq2) <- c("spp", "frequency")
region.data.freq2 <- subset(region.data.freq2, spp != "halibut")# Remove the "halibut" and "j_ak_sk" as non-target fisheries

# arrange by frequency and select the top 20
region.data.freq3 <- region.data.freq2 %>%
  dplyr::arrange(desc(frequency)) %>%
  dplyr::slice(1:15)

top_spp <- paste0(region.data.freq3$spp)
print(top_spp)

#add full names
spp_names_df <- data.frame(
  spp = top_spp,
  full_name = c("Adult Pollock","Juv. Pacific Cod","Adult Pacific Cod","Snow Crab","Adult Flathead Sole",
                "Juv. Flathead Sole","Juv. Pollock","Early Juv. Pollock","Tanner Crab","Adult Yellofin Sole",
                "Juv. Yellowfin Sole","Adult Alaska Plaice","Adult N. Rock Sole","Juv. N. Rock Sole","Juv. Alaska Skate")
)


##############################
# Fig 3 - Maximum Cohen's q
##############################

# load spearman's rank results
spearman_totals2.tmb <- readRDS("data/spearman_totals2.tmb_ar1_60_tw_m1.rds")

spearman_results <- spearman_totals2.tmb %>%
  filter(!is.na(rho)) %>%
  mutate(rho = as.numeric(rho)) %>%
  mutate(z_rho = atanh(rho)) # Fisher z-transformation 

# Filter only the dynamic models
dynamic_models <- c("annual", "timestep", "ltm")
dynamic_df <- spearman_results %>%
  filter(step %in% dynamic_models) %>%
  mutate(step = as.factor(step))

# pairwise Cohen's q
q_yearly_results <- data.frame(
  spp = character(),
  year = numeric(),
  pair = character(),
  cohens_q = numeric(),
  stringsAsFactors = FALSE
)

# function to calculate Cohen's q
cohen_q <- function(r1, r2) {
  valid_indices <- which(!is.na(r1) & !is.na(r2) & abs(r1) < 1 & abs(r2) < 1)
  
  if (length(valid_indices) == 0) {
    return(NA)
  }
  
  r1_clean <- r1[valid_indices]
  r2_clean <- r2[valid_indices]
  
  z1 <- atanh(r1_clean)
  z2 <- atanh(r2_clean)
  q <- z1 - z2
  return(q)
}

wide_df <- dynamic_df %>%
  pivot_wider(
    id_cols = c(spp, year),
    names_from = step,
    values_from = rho
  )

for (sp in top_spp) {
  sub_data <- wide_df %>% filter(spp == sp)
  combs <- combn(dynamic_models, 2)
  
  for (i in 1:ncol(combs)) {
    m1 <- combs[1, i]
    m2 <- combs[2, i]
    
    if (m1 %in% colnames(sub_data) && m2 %in% colnames(sub_data)) {
      q_values <- cohen_q(sub_data[[m1]], sub_data[[m2]])
      
      q_yearly_results <- rbind(q_yearly_results, data.frame(
        spp = sp,
        year = sub_data$year[!is.na(q_values)],
        pair = paste0(m1, "_vs_", m2),
        cohens_q = abs(q_values[!is.na(q_values)]),
        stringsAsFactors = FALSE
      ))
    }
  }
}

ttest_results <- q_yearly_results %>%
  group_by(spp, pair) %>%
  filter(n() >= 2) %>%
  summarise(
    max_q = max(cohens_q, na.rm = TRUE),
    mean_q = mean(cohens_q, na.rm = TRUE),
    t_statistic = list(t.test(cohens_q, mu = 0)),
    p_value = list(t_statistic[[1]]$p.value),
    .groups = "drop"
  ) %>%
  unnest(p_value)

# join full names
ttest_results <- ttest_results %>% dplyr::left_join(spp_names_df, by="spp")

print("Pairwise Cohen's q Results:")
print(ttest_results)

# summarize results by max_q
ttest_results2 <- ttest_results %>% dplyr::group_by(full_name) %>%
  summarise(
    max_q = max(max_q)
  )

# order high to low
ttest_results2$full_name <- fct_reorder(
  .f = ttest_results2$full_name,
  .x = ttest_results2$max_q,
  .fun = max, 
  .desc = FALSE 
)

# add column to shade above 0.1 blue
ttest_results2 <- ttest_results2 %>%
  mutate(is_significant = ifelse(max_q > 0.1, "Above Threshold", "Below Threshold"))

# upper quartile
upper_quartile <- quantile(ttest_results2$max_q, probs = 0.75, na.rm = TRUE)

# filter the species in the upper quartile
upper_species <- ttest_results2 %>%
  filter(max_q >= upper_quartile) %>%
  arrange(desc(max_q))

# add column to add upper quartile 
ttest_results2 <- ttest_results2 %>%
  mutate(upper_quart = ifelse(max_q > upper_quartile, "upper", "lower"))

# final column for colors
ttest_results2$combined_label <- paste(ttest_results2$is_significant, ttest_results2$upper_quart)

# plot
hist2 <- ggplot(ttest_results2, aes(x = full_name, y = max_q)) +
  geom_col(aes(fill = combined_label), color = NA) + 
  labs(
       x = "Species",
       y = "Cohen's q") +
  geom_hline(yintercept = 0.10, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = upper_quartile, linetype = "dashed", color = "red") +
  annotate("text", x = Inf, y = 0.1, label = "(small difference)", 
           hjust = 6.9, vjust = -0.5, color = "black", size = 3) +
  annotate("text", x = Inf, y = upper_quartile, label = "(upper quartile)", 
           hjust = 7.75, vjust = -0.5, color = "black", size = 3) +
    ylim(0, 0.4) +
  scale_fill_manual(values = c("Above Threshold lower" = "steelblue", 
                               "Above Threshold upper" = "steelblue4", 
                               "Below Threshold lower" = "gray50")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") 
hist2
ggsave(paste0("figures/Fig3_Maximum_Cohens_q.jpeg"), hist2, width=7, height=5, units="in", dpi=300)


######################################################
#set to spp above quartile for remainder of script
######################################################

diff_spp1 <- subset(ttest_results2, upper_quart=="upper")
diff_spp <- diff_spp1 %>% dplyr::left_join(spp_names_df, by="full_name")
diff_spp <- as.character(diff_spp$spp)
diff_spp


######################################################
#Fig 4 - Spearman's Rho
######################################################

#add full names
spp_names_df <- data.frame(
  spp = top_spp,
  full_name = c("Adult Pollock","Juvenile Pacific Cod","Adult Pacific Cod","Snow Crab","Adult Flathead Sole",
                "Juvenile Flathead Sole","Juvenile Pollock","Early Juvenile Pollock","Tanner Crab","Adult Yellofin Sole",
                "Juvenile Yellowfin Sole","Adult Alaska Plaice","Adult N. Rock Sole","Juvenile N. Rock Sole","Juvenile Alaska Skate")
)

spearman_totals2.tmb <- readRDS("data/spearman_totals2.tmb_ar1_60_tw_m1.rds")

#filter to spp above 0.1 rho differences
spearman_totals2.tmb <- subset(spearman_totals2.tmb, spp %in% diff_spp)
spearman_totals2.line <- spearman_totals2.tmb

#join full names
spearman_totals2.line <- spearman_totals2.line %>% dplyr::left_join(spp_names_df, by="spp")

#rename steps for plotting order
spearman_totals2.line$step <- gsub("annual$", "1: annual", spearman_totals2.line$step)
spearman_totals2.line$step <- gsub("timestep$", "2: timestep", spearman_totals2.line$step)
spearman_totals2.line$step <- gsub("ltm_static$", "4: ltm_static", spearman_totals2.line$step)
spearman_totals2.line$step <- gsub("ltm$", "3: ltm", spearman_totals2.line$step)

facet.order <- c("Adult Pacific Cod", "Juvenile Pacific Cod", 
                 "Adult Pollock", "Juvenile Pollock")

#plot in order
spearman.line.full <- ggplot(spearman_totals2.line) + 
  geom_line(aes(x = year, y = rho, color = step, linetype=step), linewidth = 0.6) + 
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5)) +
  xlab("Year") +
  labs(color = "Resolution", linetype = "Resolution") +
  scale_color_manual(labels = c("Annual", "5-yr Mean", "LTM", "Static"), 
                     values = c("red", "#E69F00", "#56B4E9", "black")) +
  scale_linetype_manual(labels = c("Annual", "5-yr Mean", "LTM", "Static"),
                        values = c("solid", "solid", "solid", "dotted")) +
  facet_wrap(~factor(full_name, levels = facet.order), scales="free", ncol=2)
spearman.line.full
ggsave(paste0("figures/Fig4_Spearmans_rho_Forecast.jpeg"), spearman.line.full, width=8, height=5, units="in", dpi=300)



######################################################
#Fig 5 
######################################################

#add full names
spp_names_df <- data.frame(
  spp = top_spp,
  full_name = c("Adult Pollock","Juvenile Pacific Cod","Adult Pacific Cod","Snow Crab","Adult Flathead Sole",
                "Juvenile Flathead Sole","Juvenile Pollock","Early Juvenile Pollock","Tanner Crab","Adult Yellofin Sole",
                "Juvenile Yellowfin Sole","Adult Alaska Plaice","Adult N. Rock Sole","Juvenile N. Rock Sole","Juvenile Alaska Skate")
)

# load spearman's rank results
spearman_totals2.tmb <- readRDS("data/spearman_totals2.tmb_ar1_60_tw_m1.rds")

#filter to spp above 0.1 rho differences
spearman_totals2.tmb <- subset(spearman_totals2.tmb, spp %in% diff_spp)

#join temp
temps <- region.data %>% dplyr::select(year, btemp)

baseline_temps <- subset(temps, year < 2009)
mean_baseline_temps <- mean(baseline_temps$btemp)

mean_annual_temps <- temps %>% dplyr::group_by(year) %>%
  summarize(
    mean_temp = mean(btemp)
  )

mean_annual_temps$temp_anom <- mean_annual_temps$mean_temp - mean_baseline_temps 
mean_annual_temps$abs_temp_anom <- abs(mean_annual_temps$temp_anom)

# join temp to spearman rhos
spearman_data <- spearman_totals2.tmb %>% dplyr::left_join(mean_annual_temps, by=c("year"))

# Ensure rho is numeric and handle NAs
spearman_data_clean <- spearman_data %>%
  filter(!is.na(rho)) %>%
  mutate(rho = as.numeric(rho)) %>%
  mutate(z_rho = atanh(rho))

# define models
dynamic_models <- c("annual", "timestep", "ltm")
static_model <- "ltm_static"

# calc the differences for each dynamic model relative to the static 
rho_differences_vs_static <- spearman_data_clean %>%
  filter(step %in% c(dynamic_models, static_model)) %>%
  pivot_wider(
    id_cols = c(year, spp), 
    names_from = step,
    values_from = z_rho,
    names_prefix = "rho_") %>%
  mutate(
    diff_annual_vs_static = rho_annual - rho_ltm_static,
    diff_timestep_vs_static = rho_timestep - rho_ltm_static,
    diff_ltm_vs_static = rho_ltm - rho_ltm_static) %>%
  
  # long format for plotting
  pivot_longer(
    cols = starts_with("diff_"),
    names_to = "comparison_type",
    values_to = "rho_difference"
  ) %>%
  
  left_join(
    spearman_data_clean %>% select(year, spp, mean_temp, temp_anom, abs_temp_anom) %>% distinct(),
    by = c("year", "spp")) %>%
  mutate(
    comparison_type = case_when(
      comparison_type == "diff_annual_vs_static" ~ "Annual vs Static",
      comparison_type == "diff_timestep_vs_static" ~ "5yr Mean vs Static",
      comparison_type == "diff_ltm_vs_static" ~ "LTM vs Static",
      TRUE ~ comparison_type
    ) %>% factor(levels = c("Annual vs Static", "5yr Mean vs Static", "LTM vs Static"))
  )


rho_differences_vs_static <- subset(rho_differences_vs_static, spp %in% diff_spp)

# join full names
rho_differences_vs_static <- rho_differences_vs_static %>% dplyr::left_join(spp_names_df, by="spp")

# set a facet order for the plot
facet.order <- c("Adult Pacific Cod", "Juvenile Pacific Cod", 
                 "Adult Pollock", "Juvenile Pollock")

# plot
plot_diff_vs_temp <- ggplot(subset(rho_differences_vs_static), aes(x = temp_anom, y = rho_difference)) +
  geom_point(aes(color = comparison_type), alpha = 0.6, size = 2) +
  geom_smooth(aes(color = comparison_type, fill = comparison_type), method = "gam", formula = y ~ s(x, k=3), se = TRUE, linewidth = 0.8, alpha = 0.25) +
  geom_hline(yintercept = 0.0, linetype = "dashed", color = "black") +
   ylim(-0.25, 0.35) +
  labs(
    x = "Mean Bottom Temperature Anomaly (°C)",
    y = "Fisher-Standardized Rho Difference (Dynamic - Static)",
    color = "Comparison"
  ) +
  scale_color_manual(values = c("Annual vs Static" = "red",
                                "5yr Mean vs Static" = "#E69F00",
                                "LTM vs Static" = "#56B4E9")) +
  scale_fill_manual(values = c("Annual vs Static" = "red",
                               "5yr Mean vs Static" = "#E69F00",
                               "LTM vs Static" = "#56B4E9"),
                    guide="none") +
  facet_wrap(~ factor(full_name, levels = facet.order), scales="free", ncol = 2, nrow=2) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10),
    legend.position = "inside",
    c(0.75, 0.25),
    legend.justification = c(0.67, 0.35),
    legend.background = element_rect(fill = "white", color = "black", linewidth=0.3),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9), 
    legend.key.size = unit(0.8, "lines") 
  )

print(plot_diff_vs_temp)

# save
ggsave(paste0("figures/Fig5_dynamic_vs_static_rho_diff.jpeg"), plot_diff_vs_temp, width=8, height=7, units="in", dpi=300)





####################################################
#  predictions
####################################################

# load prediction grid 
grid1 <- readRDS(paste0("data/grid/grid.agg10.newbathy.ebs.rds"))
grid <- grid1
max(grid$bdepth) 
max(region.data$bdepth) 
grid <- subset(grid, bdepth < max(region.data$bdepth))

#split grid into years
rm(list=ls(pattern = "grid_"))
for(j in c(1986:2019, 2021:2023)){
  assign(paste0("grid_", j), subset(grid, year==j))
}


#sdmTMB--------------------------------------------------------------------------------------------

#run all predictions

#clear files
rm(list=ls(pattern = "predictions.tmb."))
rm(list=ls(pattern = "N.model_"))

### # predictions --- (only need to run once) ----
for(spp in diff_spp) {
  for(step in c("annual", "timestep", "ltm", "ltm_static")){
    for(j in c(2009:2019, 2021:2023)){
      N.model = readRDS(paste0("species_models/tmb_ar1_60_tw_m1/tmb.fit.ebs.", step, "_", spp, ".rds"))
      pred.grid = get(paste0("grid_", j))
      pred <- stats::predict(N.model, type="response", newdata=pred.grid)
      saveRDS(pred, paste0("data/predictions/predictions.tmb.", step, "_", spp, "_", j, ".rds"))
      assign(paste0("predictions.tmb.", step, "_", spp, "_", j), pred)

    }}}

rm(list=ls(pattern = "rast1"))
rm(list=ls(pattern = "rast_"))
rm(list=ls(pattern = "rast.stack_"))

spp_rasts <- diff_spp

# main loop by species
for(spp in spp_rasts){
  for(step in c("annual", "timestep", "ltm", "ltm_static")){
    for(yr in c(2009:2019, 2021:2023)){
      # read in predictions 
      pred <- readRDS(paste0("data/predictions/predictions.tmb.", step, "_", spp, "_", yr, ".rds"))
      # convert x and y back to m
      pred$x <- as.numeric(pred$x * 100000) 
      pred$y <- as.numeric(pred$y * 100000)
      pred2 <- pred %>% dplyr::select(x, y, est)
      ras2 <- terra::rast(pred2, type="xyz", crs="epsg:3338")
      assign(paste0("rast_", spp, "_", yr, "_", step), ras2)
  }}

  # make background layer for mapping
  bnd <- ras2
  names(bnd) <- "bnd"
  bnd <- ifel(!is.na(bnd), 1, NA)
  bnd2 <- terra::as.polygons(bnd)

  # abundance plots ----------------------------------------------
  for(yr in c(2009:2019, 2021:2023)){
    
    # pull together list of all timesteps for this species and year
    rast.list <- mget(ls(pattern=paste0("rast_", spp, "_", yr)))

    #c reate raster stack for all timesteps by year
    rast_stack <- terra::rast(rast.list)

    #r ename and reorder for plots
    names(rast_stack) <- c("Annual", "Long-term Mean (LTM)", "Static", "5-year Window")
    rast_stack <- subset(rast_stack, c("Annual", "5-year Window", "Long-term Mean (LTM)", "Static"))

    # make abundance plot
    abund_plot <- ggplot() +
      tidyterra::geom_spatvector(data=bnd2, color="black", fill=NA) +
      tidyterra::geom_spatraster(data=rast_stack) +
      scale_fill_viridis_c(na.value="transparent", name="Abundance ") +
      theme_bw() + 
      facet_wrap(~factor(lyr, levels = c("Annual", "5-year Window", "Long-term Mean (LTM)", "Static")), nrow=1)

    
    # threshold plots-----------------------------------------------------------------------

    # threshold function for top 50% (the 50th percentile) 
    threshold_top50 <- function(x) {
      threshold <- quantile(x, probs = 0.5, na.rm = TRUE)
      return (x >= threshold)
    }

    # apply to each layer (makes each cell TRUE or FALSE depending on whether it meets threshold)
    rm(list=ls(pattern = "thresh_stack_"))
    for(i in 1:4){
    assign(paste0("thresh_stack_", i), app(rast_stack[[i]], threshold_top50))
    }
    
    # create stack
    stacks <- mget(ls(pattern="thresh_stack_"))
    threshold_stack <- terra::rast(stacks)
    names(threshold_stack) <- c("Annual", "5-year Window", "LTM", "Static")
    
    # classify as factor for plotting
    thresh_cat <- classify(threshold_stack, cbind(0,0))

    thresh_cat2 <- as.factor(thresh_cat)
    thresh_cat2 <- droplevels(thresh_cat2)
    levels(thresh_cat2)

    # make threshold plot
    thresh_plot <- ggplot() +
      tidyterra::geom_spatraster(data=thresh_cat2) +
      scale_fill_manual(values=c("1"= "#31688EFF"),
        labels=c("1"=""),
        name="Predicted \n50th pctl.",
        na.translate=FALSE) +
      tidyterra::geom_spatvector(data=bnd2, color="black", fill=NA) +
      guides(fill=guide_legend(na.translate=FALSE)) +
      theme_bw() +
      theme(legend.position="none") +
      facet_wrap(~factor(lyr, levels = c("Annual", "5-year Window", "LTM", "Static")), nrow=1)
  
    assign(paste0("thresh_plot_", spp, "_", yr), thresh_plot)

    spp_title2 <- ifelse(spp == "a_pcod","Adult Pacific Cod",
                              ifelse(spp == "j_pcod", "Juvenile Pacific Cod",
                                     ifelse(spp == "a_poll", "Adult Pollock", "Juvenile Pollock")))
    
    # sum across layers to get consensus
    consensus_ras <- sum(threshold_stack, na.rm = TRUE)
    
    # convert to long form for plotting
    consensus_df <- terra::as.data.frame(consensus_ras, xy=TRUE)
    consensus_df$name <- "Agreement on 50th Percentile"
    consensus_df$name2 <- spp_title2
    consensus_df_2 <- subset(consensus_df, consensus_df$sum > 0)
    
    # classify as factor for plotting
    consensus_ras2 <- consensus_ras
    consensus_ras2[consensus_ras2==0] <- NA
    consensus_ras3 <- classify(consensus_ras2, cbind(0,0))
    consensus_ras3 <- as.factor(consensus_ras3)
    consensus_ras3 <- droplevels(consensus_ras3)
    levels(consensus_ras3)
    names(consensus_ras3) <- "Agreement"
    
    # combined figure plot
    consensus_plot2 <- ggplot() +
      tidyterra::geom_spatraster(data=consensus_ras3, aes(fill=Agreement)) +
      tidyterra::geom_spatvector(data=bnd2, color="black", fill=NA) +
      scale_fill_manual(values=c(
        "0" = "transparent",
        "1" = "#000004FF",
        "2" = "#781C6DFF",
        "3" = "#ED6925FF",
        "4" = "#FCFFA4FF"),
        labels =c("0",
                  "1",
                  "2",
                  "3",
                  "4"),
        na.translate=FALSE,
        name = "Agreement \n(# Models)") +
      guides(fill=guide_legend(na.translate=FALSE)) +
      theme_bw() +
      theme(axis.title = element_blank()) +
      theme(legend.position="none") +
      facet_wrap(~lyr)
    consensus_plot2


#what percentage of total overlap? ------------------------------------------------------------------

# first calc total area (need to use rasters)
total_area <- terra::global(consensus_ras %in% c(0:4), "sum", na.rm=TRUE)

# then sum area for agreements
area_0 <- terra::global(consensus_ras == 0, "sum", na.rm=TRUE)
area_1 <- terra::global(consensus_ras == 1, "sum", na.rm=TRUE)
area_2 <- terra::global(consensus_ras == 2, "sum", na.rm=TRUE)
area_3 <- terra::global(consensus_ras == 3, "sum", na.rm=TRUE)
area_4 <- terra::global(consensus_ras == 4, "sum", na.rm=TRUE)

# finally, calc percentages of agreement
agreement_0 <- (area_0/total_area) * 100
agreement_1 <- (area_1/total_area) * 100
agreement_2 <- (area_2/total_area) * 100
agreement_3 <- (area_3/total_area) * 100
agreement_4 <- (area_4/total_area) * 100

agreement_table <- data.frame(
  num_models = 0:4,
  area = c(agreement_0[1,1], agreement_1[1,1], agreement_2[1,1], agreement_3[1,1], agreement_4[1,1]))

#now add agreement on zeros
agreement_table$spp <- spp
agreement_table$year <- yr
assign(paste0("agreement_tmb_table.", spp, ".", yr), agreement_table)

# stand-alone consensus figure plot
full_agreement <- subset(agreement_table, num_models %in% c(0,4)) %>% group_by(spp, year) %>%
  summarise(
    area = sum(area) 
  )

facet.order <- c("Adult Pacific Cod", "Juv. Pacific Cod", 
                 "Adult Pollock", "Juv. Pollock")

full_name_title <- ifelse(spp == "a_pcod","Adult Pacific Cod",
                          ifelse(spp == "j_pcod", "Juvenile Pacific Cod",
                                 ifelse(spp == "a_poll", "Adult Pollock", "Juvenile Pollock")))

#labels
new_labels <- c("Agreement" = paste0(full_name_title))

consensus_plot <- ggplot() +
  tidyterra::geom_spatraster(data=consensus_ras3, aes(fill=Agreement)) +
  tidyterra::geom_spatvector(data=bnd2, color="black", fill=NA) +
  scale_fill_manual(values=c(
    "1" = "#000004FF",
    "2" = "#781C6DFF",
    "3" = "#ED6925FF",
    "4" = "#FCFFA4FF"),
    labels =c("1",
              "2",
              "3",
              "4"),
    na.translate=FALSE,
    name = "Agreement \n(# Models)") +
  guides(fill=guide_legend(na.translate=FALSE)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(),
        ) +
  facet_wrap(~lyr, labeller=labeller(lyr=new_labels)) +
  geom_text(data=full_agreement, x=-450000, y=1400000, aes(label=paste0("Full = ", round(area,2), "%")),
            inherit.aes=FALSE, size=3)
consensus_plot
assign(paste0("consensus_plot_tmb_", spp, "_", yr), consensus_plot) # save each consensus plot

}}  

 
#pull agreement tables together
all_agreement_tables.tmb <- mget(ls(pattern="agreement_tmb_table."))
all_agreement_tables_df.tmb <- do.call(rbind, all_agreement_tables.tmb) 
all_agreement_tables_df2.tmb <- data.table::setDT(all_agreement_tables_df.tmb, keep.rownames=FALSE)[]
all_agreement_tables_df2.tmb <- subset(all_agreement_tables_df2.tmb, spp %in% spp_rasts)

all_agreement_tables_full <- all_agreement_tables_df2.tmb

all_agreement_tables_df3.tmb <- subset(all_agreement_tables_df2.tmb, num_models %in% c(0,4)) %>% group_by(spp, year) %>%
  summarise(
   area = sum(area) 
  )

all_agreement_tables_df4.tmb <- all_agreement_tables_df3.tmb %>% dplyr::group_by(spp) %>%
  reframe(
    mean_agreement = mean(area),
    sd_agreement = sd(area),
    min_agreement = min(area),
    max_agreement = max(area)
  )

all_agreement_tables_df4.tmb$range <- all_agreement_tables_df4.tmb$max_agreement - all_agreement_tables_df4.tmb$min_agreement


all_agreements <- subset(all_agreement_tables_full, num_models %in% c(0,4)) 
all_agreements2.tmb <- all_agreements %>% group_by(spp, year) %>%
  summarise(
    area = sum(area) 
  )


#############################################################
### Fig 5 - Spatial Agreement by Temperature Anomaly
#############################################################

# add full names
spp_names_df <- data.frame(
  spp = top_spp,
  full_name = c("Adult Pollock","Juvenile Pacific Cod","Adult Pacific Cod","Snow Crab","Adult Flathead Sole",
                "Juvenile Flathead Sole","Juvenile Pollock","Early Juvenile Pollock","Tanner Crab","Adult Yellofin Sole",
                "Juvenile Yellowfin Sole","Adult Alaska Plaice","Adult N. Rock Sole","Juvenile N. Rock Sole","Juvenile Alaska Skate")
)

# agreements
all_agreement_tables_temp <- all_agreements2.tmb %>% dplyr::left_join(mean_annual_temps, by=c("year"))

# join full names
all_agreement_tables_temp <- all_agreement_tables_temp %>% dplyr::left_join(spp_names_df, by="spp")

facet.order <- c("Adult Pacific Cod", "Juvenile Pacific Cod", 
                 "Adult Pollock", "Juvenile Pollock")

temp_plot_agreement <- ggplot() +
  geom_point(data=subset(all_agreement_tables_temp), aes(x=temp_anom, y=area), shape=1, size=2, color="steelblue") +
  geom_smooth(data=subset(all_agreement_tables_temp), aes(x=temp_anom, y=area), method = "gam", formula = y ~ s(x, k = 3), color="black", linewidth=0.7) +
  labs(
    x = "Mean Bottom Temperature Anomaly (°C)",
    y = "Percentage of Area Agreement"
  ) +
  theme_bw() +
  facet_wrap(~factor(full_name, levels = facet.order), scales="free") +
  theme(
    strip.text = element_text(size = 10)
  )
temp_plot_agreement
ggsave(paste0("figures/Fig6_Agreement_Temp_Anomalies.jpeg"), temp_plot_agreement, width=8, height=6, units="in", dpi=300)



#############################################################
### Fig 7 - Spatial Agreement of Pcod by Temperature Periods
#############################################################

# skill-based consensus plots
temps.ranges <- subset(mean_annual_temps, year > 2008)
temps.ranges.sum <- temps.ranges %>% 
  summarise(
    min_anom = min(temp_anom),
    min_year = year
  )
temps.ranges.sum

# years with the lowest temp_anom
lowest_anom_year <- temps.ranges %>%
  slice_min(temp_anom, n = 1, with_ties = FALSE) 

# years with the temp_anom closest to 0
closest_to_zero_anom_year <- temps.ranges %>%
  slice_min(abs(temp_anom), n = 1, with_ties = FALSE) 

# years with the largest temp_anom
largest_anom_year <- temps.ranges %>%
  slice_max(temp_anom, n = 1, with_ties = FALSE) 

# Print the results
print(lowest_anom_year) #2012 (-1.42)
print(closest_to_zero_anom_year) #2023 (-0.0001)
print(largest_anom_year) # 2019 (+2.10)

# plots

# title plot function
create_column_title_plot <- function(title_text) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = title_text,
             size = 5, fontface = "bold") +
    theme_void() +
    theme(plot.margin = margin(0,0,0,0))
}

title_col1_plot <- create_column_title_plot("Cool (2012)")
title_col2_plot <- create_column_title_plot("Average (2023)")
title_col3_plot <- create_column_title_plot("Warm (2019)")

# Print them individually:
print(title_col1_plot)
print(title_col2_plot)
print(title_col3_plot)

consensus_combined <- (title_col1_plot + title_col2_plot + title_col3_plot) / 
  (consensus_plot_tmb_a_pcod_2012 + consensus_plot_tmb_a_pcod_2023 + consensus_plot_tmb_a_pcod_2019) / 
  (consensus_plot_tmb_a_poll_2012 + consensus_plot_tmb_a_poll_2023 + consensus_plot_tmb_a_poll_2019) +
  plot_layout(
    heights = c(0.1, 1, 1),
    guides = "collect") 
consensus_combined

ggsave(paste0("figures/Fig7_temps_consensus.tmb.jpeg"), consensus_combined, width=10, height=5, units="in")





#Supplemental Figures --------------------------------------------------------------------------------------


########################################################
# Figure S1 - mesh 
########################################################
train.data = subset(region.data, year < 2009)

mesh <- sdmTMB::make_mesh(train.data, c("x", "y"),
                          fmesher_func = fmesher::fm_mesh_2d_inla,
                          cutoff = 0.6, # minimum triangle edge length
                          max.edge = c(3, 9),
                          offset=c(1.2, 9)) # inner and outer max triangle lengths
mesh$mesh$n 

#save
png("figures/Fig_S1_Mesh.png", width = 800, height = 800, res = 120)
plot(mesh)
dev.off()



##############################################################
### Figure S2 - Conditional Effects Plots
##############################################################

#annual -----------------------------------------------------------------

for(spp in top_spp) {
  for(step in c("annual")){
    #load model  
    fit <- readRDS(paste0("species_models/tmb_ar1_60_tw_m1/tmb.fit.ebs.", step, "_", spp, ".rds"))
    
    #bdepth
    bdepth_data_for_plot <- visreg(fit, xvar="bdepth_scaled", plot=FALSE)
    fit_data <- bdepth_data_for_plot$fit 
    
    # plot
    bdepth_plot <- ggplot2::ggplot(fit_data, ggplot2::aes(x = bdepth_scaled, y = visregFit)) +
      ggplot2::geom_line(color = "black") + 
      ggplot2::geom_ribbon(ggplot2::aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.5, fill = "grey") + #
      ggplot2::labs(
        x = "Bottom Depth",
        y = "f(Bottom Depth)") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "right",
                     plot.title = ggplot2::element_text(hjust = 0.5))
    
    assign(paste0("bdepth.", step, "_", spp), bdepth_plot)
    
    
    #btemp
    btemp_data_for_plot <- visreg(fit, xvar="btemp_scaled", plot=FALSE)
    
    # fitted line and CI
    fit_data <- btemp_data_for_plot$fit
    
    # plot
    btemp_plot <- ggplot2::ggplot(fit_data, ggplot2::aes(x = btemp_scaled, y = visregFit)) +
      ggplot2::geom_line(color = "black") + 
      ggplot2::geom_ribbon(ggplot2::aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.5, fill = "grey") + 
      ggplot2::labs(
        x = "Bottom Temp.",
        y = "f(Bottom Temp.)") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "right",
                     plot.title = ggplot2::element_text(hjust = 0.5))
    
    assign(paste0("btemp.", step, "_", spp), btemp_plot)
    
    #bottom currents
    current <- visreg2d(fit, x=paste0("bcurrentU_scaled"), y=paste0("bcurrentV_scaled"), plot.type='gg', nn=50)
    current <- current +
      ggplot2::geom_raster(data = current$data, ggplot2::aes(x = x, y = y, fill = z)) +
      geomtextpath::geom_textcontour(data = current$data, ggplot2::aes(x = x, y = y, z = z, label = ggplot2::after_stat(level)), color = "black") +
      ggplot2::labs(
        x = "bcurrentU", 
        y = "bcurrentV", 
        fill = "Predicted Response") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none",
                     plot.title = element_text(hjust = 0.5))
    assign(paste0("current.", step, "_", spp), current)
    
  }}


#timestep -----------------------------------------------------------------

for(spp in top_spp) {
  for(step in c("timestep")){
    #load model  
    fit <- readRDS(paste0("species_models/tmb_ar1_60_tw_m1/tmb.fit.ebs.", step, "_", spp, ".rds"))
    
    #bdepth
    bdepth_data_for_plot <- visreg(fit, xvar="bdepth_scaled", plot=FALSE)
    fit_data <- bdepth_data_for_plot$fit 
    
    # plot
    bdepth_plot <- ggplot2::ggplot(fit_data, ggplot2::aes(x = bdepth_scaled, y = visregFit)) +
      ggplot2::geom_line(color = "black") + # You can set the line color
      ggplot2::geom_ribbon(ggplot2::aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.5, fill = "grey") + 
      ggplot2::labs(
        x = "Bottom Depth",
        y = "f(Bottom Depth)") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "right",
                     plot.title = ggplot2::element_text(hjust = 0.5))
    
    assign(paste0("bdepth.", step, "_", spp), bdepth_plot)
    
    #btemp
    btemp_data_for_plot <- visreg(fit, xvar="btemp_timestep_scaled", plot=FALSE)
    fit_data <- btemp_data_for_plot$fit
    
    # plot
    btemp_plot <- ggplot2::ggplot(fit_data, ggplot2::aes(x = btemp_timestep_scaled, y = visregFit)) +
      ggplot2::geom_line(color = "black") + 
      ggplot2::geom_ribbon(ggplot2::aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.5, fill = "grey") + #
      ggplot2::labs(
        x = "Bottom Temp.",
        y = "f(Bottom Temp.)") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "right",
                     plot.title = ggplot2::element_text(hjust = 0.5))
    
    assign(paste0("btemp.", step, "_", spp), btemp_plot)
    
    #bottom currents
    current <- visreg2d(fit, x=paste0("bcurrentU_timestep_scaled"), y=paste0("bcurrentV_timestep_scaled"), plot.type='gg', nn=50)
    current <- current +
      ggplot2::geom_raster(data = current$data, ggplot2::aes(x = x, y = y, fill = z)) +
      geomtextpath::geom_textcontour(data = current$data, ggplot2::aes(x = x, y = y, z = z, label = ggplot2::after_stat(level)), color = "black") +
      ggplot2::labs(
        x = "bcurrentU", 
        y = "bcurrentV", 
        fill = "Predicted Response") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none",
                     plot.title = element_text(hjust = 0.5))
    assign(paste0("current.", step, "_", spp), current)
    
  }}


#ltm -----------------------------------------------------------------

for(spp in top_spp) {
  for(step in c("ltm")){
    #load model  
    fit <- readRDS(paste0("species_models/tmb_ar1_60_tw_m1/tmb.fit.ebs.", step, "_", spp, ".rds"))
    
    #bdepth
    bdepth_data_for_plot <- visreg(fit, xvar="bdepth_scaled", plot=FALSE)
    fit_data <- bdepth_data_for_plot$fit 
    
    # plot
    bdepth_plot <- ggplot2::ggplot(fit_data, ggplot2::aes(x = bdepth_scaled, y = visregFit)) +
      ggplot2::geom_line(color = "black") + # You can set the line color
      ggplot2::geom_ribbon(ggplot2::aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.5, fill = "grey") + 
      ggplot2::labs(
        x = "Bottom Depth",
        y = "f(Bottom Depth)") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "right",
                     plot.title = ggplot2::element_text(hjust = 0.5))
    
    assign(paste0("bdepth.", step, "_", spp), bdepth_plot)
    
    #btemp
    btemp_data_for_plot <- visreg(fit, xvar="btemp_ltm_scaled", plot=FALSE)
    fit_data <- btemp_data_for_plot$fit
    
    # plot
    btemp_plot <- ggplot2::ggplot(fit_data, ggplot2::aes(x = btemp_ltm_scaled, y = visregFit)) +
      ggplot2::geom_line(color = "black") + 
      ggplot2::geom_ribbon(ggplot2::aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.5, fill = "grey") + 
      ggplot2::labs(
        x = "Bottom Temp.",
        y = "f(Bottom Temp.)") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "right",
                     plot.title = ggplot2::element_text(hjust = 0.5))
    
    assign(paste0("btemp.", step, "_", spp), btemp_plot)
    
    #bottom currents
    current <- visreg2d(fit, x=paste0("bcurrentU_ltm_scaled"), y=paste0("bcurrentV_ltm_scaled"), plot.type='gg', nn=50)
    current <- current +
      ggplot2::geom_raster(data = current$data, ggplot2::aes(x = x, y = y, fill = z)) +
      geomtextpath::geom_textcontour(data = current$data, ggplot2::aes(x = x, y = y, z = z, label = ggplot2::after_stat(level)), color = "black") +
      ggplot2::labs(
        x = "bcurrentU", 
        y = "bcurrentV", 
        fill = "Predicted Response") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none",
                     plot.title = element_text(hjust = 0.5))
    assign(paste0("current.", step, "_", spp), current)
    
  }}


#ltm_static -----------------------------------------------------------------

for(spp in top_spp) {
  for(step in c("ltm_static")){
    #load model  
    fit <- readRDS(paste0("species_models/tmb_ar1_60_tw_m1/tmb.fit.ebs.", step, "_", spp, ".rds"))
    
    #bdepth
    bdepth_data_for_plot <- visreg(fit, xvar="bdepth_scaled", plot=FALSE)
    fit_data <- bdepth_data_for_plot$fit
    
    # plot
    bdepth_plot <- ggplot2::ggplot(fit_data, ggplot2::aes(x = bdepth_scaled, y = visregFit)) +
      ggplot2::geom_line(color = "black") + # You can set the line color
      ggplot2::geom_ribbon(ggplot2::aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.5, fill = "grey") + 
      ggplot2::labs(
        x = "Bottom Depth",
        y = "f(Bottom Depth)") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "right",
                     plot.title = ggplot2::element_text(hjust = 0.5))
    
    assign(paste0("bdepth.", step, "_", spp), bdepth_plot)
    
    #btemp
    btemp_data_for_plot <- visreg(fit, xvar="btemp_ltm_scaled", plot=FALSE)
    fit_data <- btemp_data_for_plot$fit
    
    # plot
    btemp_plot <- ggplot2::ggplot(fit_data, ggplot2::aes(x = btemp_ltm_scaled, y = visregFit)) +
      ggplot2::geom_line(color = "black") +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.5, fill = "grey") + 
      ggplot2::labs(
        x = "Bottom Temp.",
        y = "f(Bottom Temp.)") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "right",
                     plot.title = ggplot2::element_text(hjust = 0.5))
    
    assign(paste0("btemp.", step, "_", spp), btemp_plot)
    
    #bottom currents
    current <- visreg2d(fit, x=paste0("bcurrentU_ltm_scaled"), y=paste0("bcurrentV_ltm_scaled"), plot.type='gg', nn=50)
    current <- current +
      ggplot2::geom_raster(data = current$data, ggplot2::aes(x = x, y = y, fill = z)) +
      geomtextpath::geom_textcontour(data = current$data, ggplot2::aes(x = x, y = y, z = z, label = ggplot2::after_stat(level)), color = "black") +
      ggplot2::labs(
        x = "bcurrentU", 
        y = "bcurrentV", 
        fill = "Predicted Response") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none",
                     plot.title = element_text(hjust = 0.5))
    assign(paste0("current.", step, "_", spp), current)
    
  }}




#plots --------------------------------------------------------------------------------------------------

for(spp in top_spp) {
  
  #some figure names  
  title_df <- spp_names_df %>% dplyr::filter(spp_name == spp) 
  title_name <- as.character(title_df$full_name)
  
  # blank plot with only the title text
  label_annual <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, label = "Annual", size = 6) + 
    theme_void()
  
  label_5yr <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, label = "5-Year\nMean", size = 6) + 
    theme_void()
  
  label_ltm <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, label = "LTM", size = 6) + 
    theme_void()
  
  label_static <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, label = "Static", size = 6) + 
    theme_void()
  
  
  spp_plot <- ggarrange(label_annual, get(paste0("bdepth.annual_", spp)), get(paste0("btemp.annual_", spp)), get(paste0("current.annual_", spp)),
                        label_5yr, get(paste0("bdepth.timestep_", spp)), get(paste0("btemp.timestep_", spp)), get(paste0("current.timestep_", spp)),
                        label_ltm, get(paste0("bdepth.ltm_", spp)), get(paste0("btemp.ltm_", spp)), get(paste0("current.ltm_", spp)),
                        label_static, get(paste0("bdepth.ltm_static_", spp)), get(paste0("btemp.ltm_static_", spp)), get(paste0("current.ltm_static_", spp)),
                        nrow=4, ncol=4, widths = c(0.3, 1, 1, 1))
  spp_plot <- annotate_figure(spp_plot, top = ggpubr::text_grob(paste0("Conditional Effects by Temporal Resolution: ", title_name), size = 18))
  ggsave(paste0("figures/response_plots_", spp, ".jpeg"), spp_plot, width=10, height=10, units="in", dpi=300)
  
}



##############################################################
### Figure S3 - Retrospective Forecast
##############################################################

#normal lines
spearman_totals2.tmb <- readRDS("data/spearman_totals2.tmb_ar1_60_tw_m1.rds")

spearman_totals2.line <- spearman_totals2.tmb

#join full names
spearman_totals2.line <- spearman_totals2.line %>% dplyr::left_join(spp_names_df, by="spp")

#rename steps for plotting order
spearman_totals2.line$step <- gsub("annual$", "1: annual", spearman_totals2.line$step)
spearman_totals2.line$step <- gsub("timestep$", "2: timestep", spearman_totals2.line$step)
spearman_totals2.line$step <- gsub("ltm_static$", "4: ltm_static", spearman_totals2.line$step)
spearman_totals2.line$step <- gsub("ltm$", "3: ltm", spearman_totals2.line$step)

#plot in order
spearman.line.full <- ggplot(spearman_totals2.line) + 
  geom_line(aes(x = year, y = rho, color = step, linetype=step), linewidth = 0.6) + 
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5)) +
  xlab("Year") +
  ylim(0,1) +
  labs(color = "Resolution", linetype = "Resolution") +
  scale_color_manual(labels = c("Annual", "5-year Mean", "Long-term Mean (LTM)", "Static"), 
                     values = c("red", "#E69F00", "#56B4E9", "black")) +
  scale_linetype_manual(labels = c("Annual", "5-year Mean", "Long-term Mean (LTM)", "Static"),
                        values = c("solid", "solid", "solid", "dotted")) +
  facet_wrap(~factor(full_name), scales="free", ncol=3)
spearman.line.full
ggsave(paste0("figures/Fig_S3_Spearman_full.jpeg"), spearman.line.full, width=8, height=10, units="in", dpi=300)




###################################################################
### Fig S4 - Predicted Area of 50th Percentile by Temp Period
###################################################################

# --- 1. Define the Plots with Left Justification ---
create_column_title_plot <- function(title_text) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = title_text,
             size = 4,  hjust = 0.5) + 
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 10)) 
}

# title plots
title_cool <- create_column_title_plot("Cool \n(2012)")
title_average <- create_column_title_plot("Average \n(2023)")
title_warm <- create_column_title_plot("Warm \n(2019)")

# maps
map_cool <- thresh_plot_a_pcod_2012
map_average <- thresh_plot_a_pcod_2023
map_warm <- thresh_plot_a_pcod_2019

# combine title and map plots

# Row 1
row1 <- (title_cool | map_cool) + 
  plot_layout(widths = c(0.2, 1))

# Row 2
row2 <- (title_average | map_average) + 
  plot_layout(widths = c(0.2, 1))

# Row 3
row3 <- (title_warm | map_warm) + 
  plot_layout(widths = c(0.2, 1))

# stack vertically with / ---
indiv_thresh_combined <- row1 / row2 / row3

indiv_thresh_combined <- indiv_thresh_combined + 
  plot_layout(guides = "collect") 

indiv_thresh_combined
ggsave(paste0("figures/Fig_S4_temps_threshold_rasters_a_pcod.jpeg"), indiv_thresh_combined, width=8, height=5, units="in")








