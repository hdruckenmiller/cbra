# Calculate means in CBRS, counterfactuals, all coastal areas; estimation
# Date last modified: 7/7/2023
# Last modified by: Penny Liao

rm(list=ls())

library(MatchIt)
library(dplyr)
library(rgdal)
library(ggplot2)
library(PSweight) # https://cran.r-project.org/web/packages/PSweight/vignettes/vignette.pdf
library(tidyverse)
library(logr)
library(stargazer)
library(rlang)
options(warn=-1)
library(lfe)
library(estimatr)
library(formattable)
library(xtable)

date <- "07072023"
path <- "L:/Project-CBRA/"
data_path <- paste0(path, "data/")
output_path <- paste0(path, "output/")

# function for inverse hyperbolic sine transformation
ihs <- function(x) { 
  y <- log(x + sqrt(x ^ 2 + 1)) 
  return(y) 
}

# region_ids where there is a leading zero in shapefile but not match_sample
remove_leading_zero <- c("0154", "01117", "01208", "01219", "01220", "01223", 
                         "01230", "01231", "01242", "01255", "01267", "01289", "01294", 
                         "094363", "094367", "095093", "099344")

################ SETUP PROPENSITY SCORE FORMULA ################

# Load matched sample
data <- read.csv(paste0(data_path, "int/match_sample.csv"))

# Load inland counterfactual units to drop
drop <- read.csv(paste0(data_path, "int/regionalization/landuse_counterfactuals_ratio3/counterfactual_inland_drop.csv"))
names(drop) <- "region_id"
data <- data %>% filter(!region_id %in% drop$region_id)
rm(drop)

# Load barrier island/capes identifiers 
barriers <- read.csv(paste0(data_path, "int/barriers/units_barriers.csv"))
barriers$treated <- NULL
barriers$region_id[barriers$region_id %in% remove_leading_zero] <- 
  substr(barriers$region_id[barriers$region_id %in% remove_leading_zero], 2, 
         nchar(barriers$region_id[barriers$region_id %in% remove_leading_zero]))

data <- data %>% left_join(barriers, by = "region_id")

# Covariates used in match 
covariates <- c(
  # LAND COVERS
  "develpd",  "treecvr", "wetland","grssshr", 
  "croplnd", "barren", "water", "max_developed",
  
  # GEOGRAPHY 
  "elevatn", "min_dist_coast", "area_acre", 
  "region_NorthAtlantic", "region_SouthAtlantic", "region_GulfCoast",
  
  # DEVELOPMENT PRESSURES
  "pop_density", "paved_road_density", "bridge_density", 
  "pct_protected_1982", "urban_population_mean",
  
  # SOCIODEMOGRAPHICS
  "med_hhincome_p", "per_cap_income_p", "below_pov_p", 
  "ed_ba_grad_p", "employed_p", "unemployed_p",
  "female_p", "race_white_p", "race_black_p", "hispanic_p")

# Formula for PS model 
ps.any <- as.formula(paste("in_cbrs ~ ", paste(covariates, collapse=" + ")))

# Check balance 
bal.any <- SumStat(ps.formula = ps.any, data = data,
                   weight = c("overlap"))

full_sample_balance <- bal.any

summary(bal.any, metric = "PSD")
pdf("L:/Project-CBRA/output/figures/balance_plot_full_barrier.pdf", 
    width = 10, height = 8)
plot(bal.any, type = "balance", metric = "PSD")
dev.off()


# Need overlap weights for good balance
data$overlap_weights <- bal.any$ps.weights$overlap

# Saving weights
# write_csv(data[, c("region_id", "overlap_weights")], 
#           paste0(data_path, "int/est_weights/weights_fullsample_barrier.csv"))

################ Balance Check for local pushback ################

lcv <- read.csv("L:/Project-CBRA/data/raw/LCV/cbrs_controls_LCV.csv")
tax <- read.csv("L:/Project-CBRA/data/raw/PropertyTaxReliance/cbrs_control_tax.csv")

local <- merge(lcv[,c("region_id", "LCV_1982")], tax[,c("reliance", "region_id")], by = "region_id")

data <- merge(data, local, by = "region_id", all.x = T)

data_lcv <- data %>% filter(!is.na(LCV_1982))

# Full sample weights
bal.lcv <- SumStat(zname = "in_cbrs",
                      xname = c("LCV_1982", "reliance"),
                      data = data_lcv,
                      ps.estimate = data_lcv$overlap_weights,
                      trtgrp = "2",
                      weight = "overlap")

data_tax <- data %>% filter(!is.na(reliance))

# Full sample weights
bal.tax <- SumStat(zname = "in_cbrs",
                   xname = c("LCV_1982", "reliance"),
                   data = data_tax,
                   ps.estimate = data_tax$overlap_weights,
                   trtgrp = "2",
                   weight = "overlap")

local_compare <- rbind(summary(bal.lcv, metric = "PSD")$overlap[1,c(1,2,5)],
                       summary(bal.tax, metric = "PSD")$overlap[2,c(1,2,5)])
row.names(local_compare) <- c("LCV_1982", "property_tax_reliance")
local_compare

################ DIRECT EFFECT ESTIMATION ################

### Data Assembly ###

# Read land area
unit_area <- read.csv(paste0(data_path, "int/outcomes/nlcd_cbrs_counterfactuals.csv"))

# Read outcome: land values
fmv1 <- read.csv(paste0(data_path, "int/outcomes/fmv_cbrs.csv"))
fmv1 <- fmv1 %>% rename(region_id = Unit)

fmv2 <- read.csv(paste0(data_path, "int/outcomes/fmv_counterfactual.csv"))
fmv <- rbind(fmv1, fmv2)

# Read outcome: buildings
bldg1 <- read.csv(paste0(data_path, "int/outcomes/CBRS_buildings.csv"))
bldg2 <- read.csv(paste0(data_path, "int/outcomes/counterfactual_buildings.csv"))
bldg <- rbind(bldg1, bldg2)

# Read outcome: ZTRAX transactions
ztrans1 <- read.csv(paste0(data_path, "int/outcomes/ztrax_cbrs.csv"))
ztrans2 <- read.csv(paste0(data_path, "int/outcomes/ztrax_counterfactuals.csv"))
ztrans <- rbind(ztrans1, ztrans2)

# Read outcome: ZTRAX assessments
zasmt <- read.csv(paste0(data_path, "int/outcomes/zasmt_units.csv"))

# Read outcome: Census 2020 demographics
cbrs_census <- read.csv(paste0(data_path, "int/outcomes/cbrs_census2020.csv"))
counterfactual_census <- read.csv(paste0(data_path, "int/outcomes/counterfactual_census2020.csv"))
direct_census <- cbrs_census %>% rename(region_id = Unit) %>% rbind(counterfactual_census)
direct_census <- direct_census %>% dplyr::select(region_id,
                                          occupied_p_2020 = occupied_p,
                                          white_p_2020 = white_p,
                                          black_p_2020 = black_p,
                                          hispanic_p_2020 = hispanic_p)

# Fix region_ids where there is a leading zero in shapefile but not match_sample
direct_census$region_id[direct_census$region_id %in% remove_leading_zero] <- 
  substr(direct_census$region_id[direct_census$region_id %in% remove_leading_zero], 2, 
         nchar(direct_census$region_id[direct_census$region_id %in% remove_leading_zero]))

# Read outcome: ACS 2016-2020 5-year estimates
acs <- read.csv(paste0(data_path, "int/outcomes/region_acs.csv"))
names(acs) <- c("region_id", "white_acs", "black_acs", "hisp_acs", "bagrad_acs",
                "med_hh_inc_acs", "own_occup_acs", "rent_occup_acs", "occupied_acs",
                "med_rent_acs", "med_rent_inc_pct_acs")

# Fix region_ids where there is a leading zero in shapefile but not match_sample
acs$region_id[acs$region_id %in% remove_leading_zero] <- 
  substr(acs$region_id[acs$region_id %in% remove_leading_zero], 2, 
         nchar(acs$region_id[acs$region_id %in% remove_leading_zero]))

# Read 1982 building counts
bldg1982 <- read.csv("L:/Project-CBRA/data/int/outcomes/zasmt_pre1982_buildings.csv")

# Fix region_ids where there is a leading zero in shapefile but not match_sample
bldg1982$region_id[bldg1982$region_id %in% remove_leading_zero] <- 
  substr(bldg1982$region_id[bldg1982$region_id %in% remove_leading_zero], 2, 
         nchar(bldg1982$region_id[bldg1982$region_id %in% remove_leading_zero]))

# Read LCMAP outcome
lcmap <- read.csv("L:/Project-CBRA/data/int/outcomes/LCMAP_landcover_2018.csv")

# Fix region_ids where there is a leading zero in shapefile but not match_sample
lcmap$region_id[lcmap$region_id %in% remove_leading_zero] <- 
  substr(lcmap$region_id[lcmap$region_id %in% remove_leading_zero], 2, 
         nchar(lcmap$region_id[lcmap$region_id %in% remove_leading_zero]))

# Read SFHA building density outcomes
sfha_bldgs <- read.csv("L:/Project-CBRA/data/int/sfha_dev_density/sfha_dev_density.csv")

# Fix region_ids where there is a leading zero in shapefile but not match_sample
sfha_bldgs$region_id[sfha_bldgs$region_id %in% remove_leading_zero] <- 
  substr(sfha_bldgs$region_id[sfha_bldgs$region_id %in% remove_leading_zero], 2, 
         nchar(sfha_bldgs$region_id[sfha_bldgs$region_id %in% remove_leading_zero]))

# Assemble estimation dataset
est_data <- data %>% 
  left_join(fmv, by = "region_id") %>%
  left_join(bldg, by = "region_id") %>%
  left_join(ztrans, by = "region_id") %>%
  left_join(zasmt, by = "region_id") %>%
  left_join(unit_area[,c("region_id", "land_acres")], by = "region_id") %>%
  left_join(direct_census, by = "region_id") %>%
  left_join(acs, by = "region_id") %>%
  left_join(bldg1982, by = "region_id") %>%
  left_join(lcmap, by = "region_id") %>%
  left_join(sfha_bldgs, by = "region_id")

rm(fmv1, fmv2, fmv, bldg1, bldg2, bldg, ztrans1, ztrans2, ztrans, zasmt, unit_area, 
   direct_census, cbrs_census, counterfactual_census, bldg1982)

est_data <- est_data %>% dplyr::mutate(nsales = tidyr::replace_na(nsales, 0),
                                nsales_acre = nsales/land_acres,
                                TotalAssessedValue_acre = TotalAssessedValue_SUM/land_acres,
                                LandAssessedValue_acre = LandAssessedValue_SUM/land_acres,
                                ImprovementAssessedValue_acre = ImprovementAssessedValue_SUM/land_acres,
                                log_TotalAssessedValue_acre = ihs(TotalAssessedValue_acre),
                                log_LandAssessedValue_acre = ihs(LandAssessedValue_acre),
                                log_ImprovementAssessedValue_acre = ihs(ImprovementAssessedValue_acre),
                                TaxAmount_acre = TaxAmount_SUM/land_acres,
                                PropertyCount_acre = PropertyCount_SUM/land_acres,
                                pre1982_buildings_acre = pre1982_buildings/land_acres)

est_data <- est_data %>% mutate(region = case_when(region_NorthAtlantic == 1 ~ 1,
                                                   region_SouthAtlantic == 1 ~ 2,
                                                   region_GulfCoast == 1 ~ 3))

est_data$pre1982_buildings_acre <- tidyr::replace_na(est_data$pre1982_buildings_acre, 0)

### Check balance in ZTRAX sample ###

ztrax_data <- est_data %>% filter(!is.na(TotalAssessedValue_SUM))

# Full sample weights
# bal.ztrax <- SumStat(zname = "in_cbrs",
#                      xname = covariates,
#                      data = ztrax_data,
#                      ps.estimate = ztrax_data$overlap_weights,
#                      trtgrp = "2",
#                      weight = "overlap")
# summary(bal.ztrax, metric = "PSD")
# plot(bal.ztrax, type = "balance", metric = "PSD")


# Sample-specific weights
# Formula for PS model 
ps.any <- as.formula(paste("in_cbrs ~ ", paste(covariates, collapse=" + ")))

bal.any <- SumStat(ps.formula = ps.any, 
                   data = ztrax_data,
                   weight = "overlap")
summary(bal.any, metric = "PSD")
# pdf("L:/Project-CBRA/output/figures/balance_plot_ztrax_barrier.pdf", 
#     width = 10, height = 8)
# plot(bal.any, type = "balance", metric = "PSD")
# dev.off()

# saving weights
ztrax_data$overlap_weights <- bal.any$ps.weights$overlap

# write_csv(ztrax_data[, c("region_id", "overlap_weights")], 
#           paste0(data_path, "int/est_weights/weights_main_zsample_barrier.csv"))



################ Balance Check for 1982 Building Counts ################

# data_bldg82 <- est_data %>% filter(!is.na(pre1982_buildings_acre))
# 
# # Full sample weights
# bal.bldg82 <- SumStat(zname = "in_cbrs",
#                       xname = c("pre1982_buildings_acre", "YearBuilt_missing"),
#                       data = data_bldg82,
#                       ps.estimate = data_bldg82$overlap_weights,
#                       trtgrp = "2",
#                       weight = "overlap")
# summary(bal.bldg82, metric = "PSD")
# plot(bal.bldg82, type = "balance", metric = "PSD")



### Estimation using regressions - robust SE ###


# Function to run regressions 
est_felm <- function(depvar) {
  
  log_print(paste0("---------------- ", depvar, " ----------------"))
  
  # list of outcomes to use full sample weight
  depvar_full <- c("n_buildings_acre", "nsales_acre",
                   "occupied_p_2020", "white_p_2020", "black_p_2020", "hispanic_p_2020",
                   "white_acs", "black_acs", "hisp_acs", "bagrad_acs", "med_hh_inc_acs", 
                   "own_occup_acs", "rent_occup_acs", "occupied_acs",
                   "med_rent_acs", "med_rent_inc_pct_acs", "pre1982_buildings_acre",
                   "developed_2018", "cropland_2018", "treecover_2018", "wetland_2018", "barren_2018", "grassshrub_2018")
  
  # select estimation sample
  if (depvar %in% depvar_full) {
    temp_data <- est_data %>% filter(!is.na(!!sym(depvar)))
  } else {
    temp_data <- ztrax_data %>% filter(!is.na(!!sym(depvar)))
  }
  
  # regression formula   
  form <- as.formula(sprintf(paste0(depvar, " ~ in_cbrs + barrier")))
  
  # weighted least square
  wls <- lm_robust(formula = form,
                    data = temp_data,
                    weights = temp_data$overlap_weights)
  
  log_print(summary(wls))
  
  temp <- wls %>% tidy(conf.int = T) %>% filter(term == "in_cbrs")

  temp <- temp %>% mutate(outcome_mean = mean(temp_data[, depvar], na.rm = T))
  return(temp)
}



# Open log
log_open(paste0(output_path, "est_direct_", date, "_robust_se_barrier_wt.log"))

log_print("---------------- Summary stats ----------------")
log_print(stargazer(est_data, type = "text"))


# Apply to list of outcomes

outcomes <- c("n_buildings_acre", "developed_2018", "TotalAssessedValue_acre", 
              "mean_price", "LotSizeAcres_MEAN", "YearBuilt_MEAN", "sqfeet_MEAN", "TotalBedrooms_MEAN", 
              "white_acs", "black_acs", "med_hh_inc_acs", "hisp_acs", "bagrad_acs",  
              "own_occup_acs", "rent_occup_acs", "occupied_acs", "med_rent_acs", "med_rent_inc_pct_acs")

# # All possible outcomes
# outcomes <- c("n_buildings_acre", "mean_price", "nsales_acre",
#               "TotalAssessedValue_acre", "LandAssessedValue_acre", "ImprovementAssessedValue_acre", 
#               "LotSizeAcres_MEAN", "YearBuilt_MEAN", "sqfeet_MEAN", "TotalBedrooms_MEAN", 
#               "white_acs", "black_acs", "hisp_acs", "bagrad_acs", "med_hh_inc_acs", 
#               "own_occup_acs", "rent_occup_acs", "occupied_acs",
#               "med_rent_acs", "med_rent_inc_pct_acs", "pre1982_buildings_acre",
#               "developed_2018", "cropland_2018", "treecover_2018", "wetland_2018", "barren_2018", "grassshrub_2018",
#               "no_units_acre_90", "med_rent_p_90", "med_value_p_90", "no_units_acre_20", "med_rent_p_20", "med_value_p_20")

est_summ <- lapply(outcomes, est_felm)


est_summ <- bind_rows(est_summ)
est_summ <- est_summ[,c(9:10, 2:8)]

# Estimation table
rownames(est_summ) <- est_summ$outcome
stargazer(est_summ[, c(2:4, 6, 9)],
          type = "text",
          #type = "latex",
          #out = paste0(output_path, "estimation/est_direct.tex"),
          summary = FALSE)

obs <- est_summ$df
est_summ <- est_summ[,c(1:7)] %>% 
  mutate(across(where(is.numeric), ~ comma(., digits = 3))) %>%
  cbind(obs)

rownames(est_summ) <- NULL


log_print("---------------- Estimates ----------------")
log_print(est_summ)

# Close log
log_close()





### Estimation using regressions - bootstrap SE ###

# Function to run regressions 
est_felm <- function(depvar) {
  
  log_print(paste0("---------------- ", depvar, " ----------------"))
  
  # list of outcomes to use full sample weight
  depvar_full <- c("n_buildings_acre", "nsales_acre",
                   "occupied_p_2020", "white_p_2020", "black_p_2020", "hispanic_p_2020",
                   "white_acs", "black_acs", "hisp_acs", "bagrad_acs", "med_hh_inc_acs", 
                   "own_occup_acs", "rent_occup_acs", "occupied_acs",
                   "med_rent_acs", "med_rent_inc_pct_acs", "pre1982_buildings_acre",
                   "developed_2018", "cropland_2018", "treecover_2018", "wetland_2018", "barren_2018", "grassshrub_2018")
  
  
  # select estimation sample
  if (depvar %in% depvar_full) {
    temp_data <- est_data %>% filter(!is.na(!!sym(depvar)))
  } else {
    temp_data <- ztrax_data %>% filter(!is.na(!!sym(depvar)))
  }
  
  # regression formula   
  form <- as.formula(sprintf(paste0(depvar, " ~ in_cbrs + barrier"))) 
  
  # weighted least square
  wls <- felm(formula = form,
              data = temp_data,
              weight = temp_data$overlap_weights)
  
  # bootstrap SE
  n <- 1000
  boot_coefs <- rep(NA,n)
  for(i in 1:n){
    samp <- sample_n(temp_data, nrow(temp_data), replace = T)
    boot_wls <- felm(formula = form,
                     data = samp,
                     weight = samp$overlap_weights)
    boot_coefs[i] <- summary(boot_wls)$coefficients["in_cbrs",1]
  }
  se <- sqrt(sum((boot_coefs - mean(boot_coefs))^2)*(1/(n-1)))
  
  # log_print(summary(wls))
  
  temp <- wls %>% tidy(conf.int = T) %>% filter(term == "in_cbrs")
  temp$std.error <- se
  temp$statistic <- temp$estimate/temp$std.error
  temp$p.value <- 2*pt(abs(temp$statistic), summary(wls)$df[1], lower.tail = F)
  
  temp <- temp %>% select(estimate, std.error, statistic, p.value) %>%
    mutate(DepVar = depvar, DepVarMean = mean(temp_data[, depvar], na.rm = T))
  return(temp)
}



# Open log
log_open(paste0(output_path, "est_direct_", date, "_bootstrap_se.log"))

log_print("---------------- Summary stats ----------------")
log_print(stargazer(est_data, type = "text"))


# Apply to list of outcomes
outcomes <- c("n_buildings_acre", "developed_2018", "TotalAssessedValue_acre", 
              "mean_price", "LotSizeAcres_MEAN", "YearBuilt_MEAN", "sqfeet_MEAN", "TotalBedrooms_MEAN", 
              "white_acs", "black_acs", "med_hh_inc_acs", "hisp_acs", "bagrad_acs",  
              "own_occup_acs", "rent_occup_acs", "occupied_acs", "med_rent_acs", "med_rent_inc_pct_acs")

est_summ <- lapply(outcomes, est_felm)
est_summ <- bind_rows(est_summ)
est_summ <- est_summ %>% select(DepVar, DepVarMean, estimate, std.error, statistic, p.value)
est_summ <- est_summ %>% mutate(across(where(is.numeric), ~ comma(., digits = 3)))

log_print("---------------- Estimates ----------------")
log_print(est_summ)

# Close log
log_close()


# Estimation table
row.names(est_summ) <- est_summ$DepVar
est_summ$estimate <- format(round(est_summ$estimate,3),nsmall=3)
est_summ$std.error <- format(round(est_summ$std.error,3),nsmall=3)
est_summ$p.value <- format(round(est_summ$p.value,3),nsmall=3)
stargazer(est_summ[, c(1,3:4,6)],
          type = "latex",
          out = paste0(output_path, "estimation/est_direct_bootstrap_se.tex"),
          digits = 3,
          summary = FALSE)




####  Summary stats table for direct effect sample

data_summ <- est_data %>% group_by(in_cbrs) %>%
  summarize(n_buildings_acre = weighted.mean(n_buildings_acre, overlap_weights, na.rm = T),
            nsales_acre = weighted.mean(nsales_acre, overlap_weights, na.rm = T),
            white_acs = weighted.mean(white_acs, overlap_weights, na.rm = T),
            black_acs = weighted.mean(black_acs, overlap_weights, na.rm = T),
            hisp_acs = weighted.mean(hisp_acs, overlap_weights, na.rm = T),
            bagrad_acs = weighted.mean(bagrad_acs, overlap_weights, na.rm = T),
            med_hh_inc_acs = weighted.mean(med_hh_inc_acs, overlap_weights, na.rm = T),
            own_occup_acs = weighted.mean(own_occup_acs, overlap_weights, na.rm = T),
            rent_occup_acs = weighted.mean(rent_occup_acs, overlap_weights, na.rm = T),
            occupied_acs = weighted.mean(occupied_acs, overlap_weights, na.rm = T),
            med_rent_acs = weighted.mean(med_rent_acs, overlap_weights, na.rm = T),
            med_rent_inc_pct_acs = weighted.mean(med_rent_inc_pct_acs, overlap_weights, na.rm = T),
            developed_2018 = weighted.mean(developed_2018, overlap_weights, na.rm = T))
data_summ <- t(data_summ)
          

temp <- ztrax_data %>% group_by(in_cbrs) %>%
  summarize(mean_price = weighted.mean(mean_price, overlap_weights, na.rm = T),
            TotalAssessedValue_acre = weighted.mean(TotalAssessedValue_acre, overlap_weights, na.rm = T),
            LandAssessedValue_acre = weighted.mean(LandAssessedValue_acre, overlap_weights, na.rm = T),
            ImprovementAssessedValue_acre = weighted.mean(ImprovementAssessedValue_acre, overlap_weights, na.rm = T),
            LotSizeAcres_MEAN = weighted.mean(LotSizeAcres_MEAN, overlap_weights, na.rm = T),
            YearBuilt_MEAN = weighted.mean(YearBuilt_MEAN, overlap_weights, na.rm = T),
            sqfeet_MEAN = weighted.mean(sqfeet_MEAN, overlap_weights, na.rm = T),
            TotalBedrooms_MEAN = weighted.mean(TotalBedrooms_MEAN, na.rm = T),
            pre1982_buildings_acre = weighted.mean(pre1982_buildings_acre, na.rm = T))

temp <- t(temp)

data_summ <- rbind(data_summ[1:3,], temp[2:10,], data_summ[4:14,])

stargazer(data_summ, type = "latex", 
          out = "../output/tables/summstat_direct.tex")








################ SPILLOVER EFFECT ESTIMATION ################

### Data Assembly ##

# Specify spillover type
stnd_buffer <- TRUE

if (stnd_buffer == FALSE) {
  
  # Read land area
  sp_area <- read.csv(paste0(data_path, "int/outcomes/nlcd_spillovers.csv"))
  sp_area <- sp_area %>% arrange(region_id, -barrier) %>% 
    group_by(region_id) %>% filter(row_number() == 1)
  
  # Read outcome: land values
  fmv <- read.csv(paste0(data_path, "int/outcomes/fmv_spillovers_extended.csv"))
  
  # Read outcome: building footprints 
  buildings <- read.csv(paste0(data_path, "int/outcomes/spillover_buildings.csv"))
  buildings <- buildings %>% arrange(region_id, -barrier) %>%
    group_by(region_id) %>% filter(row_number() == 1) %>%
    select(-barrier)
  
  # Read outcome: ZTRAX transactions
  
  # Read outcome: ZTRAX assessments
  zasmt <- read.csv(paste0(data_path, "int/outcomes/zasmt_spillovers.csv"))
  
  zasmt <- zasmt %>% arrange(region_id, -barrier) %>%
    group_by(region_id) %>% filter(row_number() == 1) %>%
    select(-barrier)
  
} else {
  
  # Read land area
  sp_area <- read.csv(paste0(data_path, "int/nlcd/spill_nlcd.csv"))
  sp_area <- sp_area %>% mutate(land_acres = (total_nlcd_area_m2 - water_m2)* 0.00024711)
  sp_area <- sp_area %>% group_by(regin_d) %>% summarize(land_acres = sum(land_acres, na.rm = T))
  names(sp_area)[1] <- "region_id"
  
  # Read outcome: land values
  fmv <- read.csv(paste0(data_path, "int/outcomes/fmv_spillovers_2kmBuffer.csv"))
  
  # Read outcome: building footprints
  buildings <- read.csv(paste0(data_path, "int/outcomes/spillover_2km_buildings.csv"))
  buildings <- buildings %>% rename(n_buildings = no_buildings)
  
  # Read outcome: ZTRAX transactions
  ztrans <- read.csv(paste0(data_path, "int/outcomes/spillover_2km_ztrans_unit_with_logs.csv"))
  
  # Read outcome: ZTRAX assessments
  zasmt <- read.csv(paste0(data_path, "int/outcomes/spillover_2km_zasmt_unit_with_logs.csv"))
  
  # Read 2020 Census demographics
  spillover_census <- read_csv(paste0(data_path, "int/outcomes/spillover_2km_census2020.csv"))
  
  spillover_census <- spillover_census %>% dplyr::select(region_id,
                                                  occupied_p_2020 = occupied_p,
                                                  white_p_2020 = white_p,
                                                  black_p_2020 = black_p,
                                                  hispanic_p_2020 = hispanic_p)
  
  # Read NFIP outcomes
  nfip <- read_csv(paste0(data_path, "int/nfip_outcomes/nfip_spillover_2km.csv"))
  
  nfip <- nfip %>% dplyr::select(region_id, tot_claims_per_acre, tot_claims_per_1000_coverage, tot_bldg_per_acre)
  
  # Read outcome: ACS demographics
  acs <- read_csv(paste0(data_path, "int/outcomes/spillover_2km_aggregate_acs.csv"))
  
}

# Assemble estimation dataset
est_data <- fmv %>%
  left_join(buildings, by = "region_id") %>%
  left_join(ztrans, by = "region_id") %>%
  left_join(zasmt, by = "region_id") %>%
  left_join(sp_area, by = "region_id") %>%
  left_join(spillover_census, by = "region_id") %>%
  left_join(nfip, by = "region_id") %>%
  left_join(acs, by = "region_id")
  

# Fix region_ids where there is a leading zero in shapefile but not match_sample
est_data$region_id[est_data$region_id %in% remove_leading_zero] <- 
  substr(est_data$region_id[est_data$region_id %in% remove_leading_zero], 2, 
         nchar(est_data$region_id[est_data$region_id %in% remove_leading_zero]))

est_data <- left_join(data, est_data, by = "region_id")

est_data <- est_data %>% mutate(n_buildings_acre = n_buildings/land_acres,
                                sales_count_per_acre = replace_na(sales_count_per_acre, 0),
                                log_total_assessed_value_per_acre = ihs(total_assessed_value_per_acre),
                                log_land_assessed_value_per_acre = ihs(total_land_assessed_value_per_acre),
                                log_improvement_assessed_value_per_acre = ihs(total_improvement_assessed_value_per_acre))

rm(fmv, zasmt, sp_area, buildings, ztrans, spillover_census, acs, nfip)


### Check balance in ZTRAX sample ###

ztrax_data <- est_data %>% filter(!is.na(total_assessed_value))

# # Full sample weights
# bal.ztrax <- SumStat(zname = "in_cbrs",
#                      xname = covariates,
#                      data = ztrax_data,
#                      ps.estimate = ztrax_data$overlap_weights,
#                      trtgrp = "2",
#                      weight = "overlap")
# summary(bal.ztrax, metric = "PSD")
# plot(bal.ztrax, type = "balance", metric = "PSD")


# Sample-specific weights
bal.any <- SumStat(ps.formula = ps.any, 
                   data = ztrax_data,
                   weight = "overlap")
summary(bal.any, metric = "PSD")
plot(bal.any, type = "balance", metric = "PSD")

ztrax_sample_bal <- bal.any

# saving weights
ztrax_data$overlap_weights <- bal.any$ps.weights$overlap

#write_csv(ztrax_data[, c("region_id", "overlap_weights")], 
#          paste0(data_path, "int/est_weights/weights_stndbuffer_zsample.csv"))


### Estimation using regressions ###

# Function to run regressions (same output as above)
est_felm <- function(depvar) {
  
  log_print(paste0("---------------- ", depvar, " ----------------"))
  
  # filter for non-missing outcome
  if (depvar %in% c("fmv_lv", "n_buildings_acre", 
                    "occupied_p_2020", "white_p_2020", "black_p_2020", "hispanic_p_2020")) {
    temp_data <- est_data %>% filter(!is.na(!!sym(depvar)))
  } else {
    temp_data <- ztrax_data %>% filter(!is.na(!!sym(depvar)))
  }
  
  # regression formula   
  form <- as.formula(sprintf(paste0(depvar, " ~ in_cbrs + barrier")))
  
  # weighted least square
  wls <- lm_robust(formula = form,
                   data = temp_data,
                   weights = temp_data$overlap_weights)
  
  log_print(summary(wls))
  
  temp <- wls %>% tidy(conf.int = T) %>% filter(term == "in_cbrs")
  
  temp <- temp %>% mutate(outcome_mean = mean(temp_data[, depvar], na.rm = T))
  return(temp)
}


# Open log
if (stnd_buffer == FALSE) {
  log_open(paste0(output_path, "est_spillovers_extended_", date, "_reg.log"))
} else {
  log_open(paste0(output_path, "est_spillovers_2km_", date, "_robust_se.log"))
}

log_print("---------------- Summary stats ----------------")
log_print(stargazer(est_data, type = "text"))


# Apply to list of outcomes
outcomes <- c(
  "fmv_lv", "n_buildings_acre", "avg_sales_price_amount", "sales_count_per_acre",
  "total_assessed_value_per_acre", "total_improvement_assessed_value_per_acre", "total_land_assessed_value_per_acre",  
  "avg_lot_size", "avg_year_built", "avg_sqfeet", "avg_bedrooms",
  "occupied_p_2020", "white_p_2020", "black_p_2020", "hispanic_p_2020"
  )


est_summ <- lapply(outcomes, est_felm)
est_summ <- bind_rows(est_summ)
est_summ <- est_summ[,c(9:10, 2:8)]

# Estimation table
rownames(est_summ) <- est_summ$outcome
stargazer(est_summ[, c(2:4, 6, 9)],
          type = "latex",
          out = paste0(output_path, "estimation/est_spillover_2km.tex"),
          summary = FALSE)

obs <- est_summ$df
est_summ <- est_summ[,c(1:7)] %>% 
  mutate(across(where(is.numeric), ~ comma(., digits = 3))) %>%
  cbind(obs)
rownames(est_summ) <- NULL


log_print("---------------- Estimates ----------------")
log_print(est_summ)

# Close log
log_close()


####  Summary stats table for direct effect sample

data_summ <- est_data %>% group_by(in_cbrs) %>%
  summarize(n_buildings_acre = weighted.mean(n_buildings_acre, overlap_weights, na.rm = T),
            sales_count_per_acre = weighted.mean(sales_count_per_acre, overlap_weights, na.rm = T),
            tot_claims_per_acre = weighted.mean(tot_claims_per_acre, overlap_weights, na.rm = T),
            tot_claims_per_1000_coverage = weighted.mean(tot_claims_per_1000_coverage, overlap_weights, na.rm = T),
            tot_bldg_per_acre = weighted.mean(tot_bldg_per_acre, overlap_weights, na.rm = T),
            white = weighted.mean(white, overlap_weights, na.rm = T),
            black = weighted.mean(black, overlap_weights, na.rm = T),
            hispanic = weighted.mean(hispanic, overlap_weights, na.rm = T),
            bagrad = weighted.mean(bagrad, overlap_weights, na.rm = T),
            median_hh_income = weighted.mean(median_hh_income, overlap_weights, na.rm = T),
            own_occupied = weighted.mean(own_occupied, overlap_weights, na.rm = T),
            rent_occupied = weighted.mean(rent_occupied, overlap_weights, na.rm = T),
            occupied = weighted.mean(occupied, overlap_weights, na.rm = T),
            median_rent = weighted.mean(median_rent, overlap_weights, na.rm = T),
            median_rent_inc_pct = weighted.mean(median_rent_inc_pct, overlap_weights, na.rm = T))
data_summ <- t(data_summ)


temp <- ztrax_data %>% group_by(in_cbrs) %>%
  summarize(avg_sales_price_amount = weighted.mean(avg_sales_price_amount, overlap_weights, na.rm = T),
            total_assessed_value_per_acre = weighted.mean(total_assessed_value_per_acre, overlap_weights, na.rm = T),
            total_land_assessed_value_per_acre = weighted.mean(total_land_assessed_value_per_acre, overlap_weights, na.rm = T),
            total_improvement_assessed_value_per_acre = weighted.mean(total_improvement_assessed_value_per_acre, overlap_weights, na.rm = T),
            avg_lot_size = weighted.mean(avg_lot_size, overlap_weights, na.rm = T),
            avg_year_built = weighted.mean(avg_year_built, overlap_weights, na.rm = T),
            avg_sqfeet = weighted.mean(avg_sqfeet, overlap_weights, na.rm = T),
            avg_bedrooms = weighted.mean(avg_bedrooms, na.rm = T))

temp <- t(temp)

data_summ <- rbind(data_summ[1:3,], temp[2:9,], data_summ[4:16,])

stargazer(data_summ, type = "latex", 
          out = "../output/tables/summstat_spillover.tex")








################## BALANCE TABLES AND FIGURES #############

bal.any <- full_sample_balance 

smd_overlap <- summary(bal.any, metric = "PSD")$overlap[,"SMD"]
smd_unweight <- summary(bal.any, metric = "PSD")$unweighted[,"SMD"]

order <- c("barren", "croplnd", "develpd", "max_developed", 
           "grssshr", "treecvr", "water", "wetland", 
           "elevatn", "min_dist_coast", "area_acre",
           "region_NorthAtlantic", "region_SouthAtlantic", "region_GulfCoast",
           "pop_density", "paved_road_density", "bridge_density",
           "urban_population_mean", "pct_protected_1982", 
           "med_hhincome_p", "per_cap_income_p", "below_pov_p",
           "ed_ba_grad_p", "employed_p", "unemployed_p", 
           "female_p", "race_white_p", "race_black_p", "hispanic_p")
full_names <- c("Barren", "Cropland", "Developed", 
                "Grass/shrub", "Tree cover", "Water", "Wetland", 
                "Most densely developed", 
                "Elevation", "Distance to coast", "Area of unit",
                "Region: North", "Region: South", "Region: Gulf",
                "Population density", "Raved road density", "Bridge density",
                "Urban population", "Protected area", 
                "Median household income", "Per capita income", "Below poverty line",
                "College educated", "Employed", "Unemployed", 
                "Female", "White", "Black", "Hispanic")

smd_overlap <- smd_overlap[order]
smd_overlap <- as.data.frame(smd_overlap)
colnames(smd_overlap) <- "SMD"
rownames(smd_overlap) <- full_names
smd_overlap$index <- nrow(smd_overlap):1

smd_unweight <- smd_unweight[order]
smd_unweight <- as.data.frame(smd_unweight)
colnames(smd_unweight) <- "SMD"
rownames(smd_unweight) <- full_names
smd_unweight$index <- nrow(smd_unweight):1

pdf("L:Project-CBRA/output/figures/balance_plots_clean1.pdf",
    width = 8, height = 11)
par(mfrow=c(1,2))
plot(main = "Full sample", 
     smd_overlap$SMD, 
     smd_overlap$index, xlim = c(0, 0.55), 
     pch = 1, col = "blue", 
     yaxt="n", ylab = NA, xlab = "SMD")
axis(2, at=nrow(smd_overlap):1, row.names(smd_overlap), las = 1)
abline(v = 0.1, lty = 2)
points(smd_unweight, pch = 2, col = "darkgreen")



bal.any <- ztrax_sample_bal

smd_overlap <- summary(bal.any, metric = "PSD")$overlap[,"SMD"]
smd_unweight <- summary(bal.any, metric = "PSD")$unweighted[,"SMD"]

smd_overlap <- smd_overlap[order]
smd_overlap <- as.data.frame(smd_overlap)
colnames(smd_overlap) <- "SMD"
rownames(smd_overlap) <- full_names
smd_overlap$index <- nrow(smd_overlap):1

smd_unweight <- smd_unweight[order]
smd_unweight <- as.data.frame(smd_unweight)
colnames(smd_unweight) <- "SMD"
rownames(smd_unweight) <- full_names
smd_unweight$index <- nrow(smd_unweight):1


plot(main = "ZTRAX sample", 
     smd_overlap$SMD, 
     smd_overlap$index, xlim = c(0, 0.55), 
     pch = 1, col = "blue", 
     yaxt="n", ylab = NA, xlab = "SMD")
abline(v = 0.1, lty = 2)
points(smd_unweight, pch = 2, col = "darkgreen")
dev.off()



####  Balance table 

# Mean characteristics: show unweighted means
balance_table <- as.data.frame(summary(bal.any, metric = "PSD")$unweighted)
# Add sample size 
sample_size <- as.data.frame(matrix(NA,2,5))
row.names(sample_size) <- c("Observations", "Effective Sample Size")
sample_size[1,] <- c(summary(bal.any, metric = "PSD")$effective.sample.size[,"unweighted"], NA, NA, NA)
sample_size[2,] <- c(summary(bal.any, metric = "PSD")$effective.sample.size[,"overlap"], NA, NA, NA)
colnames(sample_size) <- colnames(balance_table)
balance_table <- rbind(balance_table, sample_size)
# Reorder columns 
balance_table <- balance_table[,c(2,1)]
colnames(balance_table) <- c("CBRS_unweighted_mean", "control_unweighted_mean")
# Add SMD (population SD using weighted data)
balance_table$SMD <- c(summary(bal.any, metric = "PSD")$overlap[,"SMD"], NA, NA)
# Scale units
balance_table[c(1:8,12:14,18,22:29), 1:2] <- balance_table[c(1:8,12:14,18,22:29), 1:2]*100 # to %
balance_table["paved_road_density", 1:2] <- balance_table["paved_road_density", 1:2]*1000 # to km
balance_table["urban_population_mean", 1:2] <- balance_table["urban_population_mean", 1:2]/1000000 # to millions
balance_table["med_hhincome_p", 1:2] <- balance_table["med_hhincome_p", 1:2]/1000 # to thousands
balance_table <- balance_table[-c(8,11:14,21)]
xtable(balance_table, digits = 3)


####  Check balance in barrier island subgroups

# barrier island sample
df_barrier <- data %>% filter(barrier == 1)

bal.barrier <- SumStat(zname = "in_cbrs",
                       xname = covariates,
                       data = df_barrier,
                       ps.estimate = df_barrier$overlap_weights,
                       trtgrp = "2",
                       weight = "overlap")
summary(bal.barrier, metric = "PSD")

pdf("L:/Project-CBRA/output/figures/balance_barrier_nowt.pdf", 
    width = 10, height = 8)
plot(bal.barrier, type = "balance", metric = "PSD")
dev.off()

# mainland sample
df_mainland <- data %>% filter(barrier == 0)

bal.mainland <- SumStat(zname = "in_cbrs",
                        xname = covariates,
                        data = df_mainland,
                        ps.estimate = df_mainland$overlap_weights,
                        trtgrp = "2",
                        weight = "overlap")
summary(bal.mainland, metric = "PSD")

pdf("L:/Project-CBRA/output/figures/balance_mainland_nowt.pdf", 
    width = 10, height = 8)
plot(bal.mainland, type = "balance", metric = "PSD")
dev.off()

# weights comparison
df_mainland %>%
  ggplot() +
  geom_histogram(aes(x = overlap_weights, y = after_stat(count / sum(count))), 
                 fill = "#69b3a2", alpha = 0.7) + 
  geom_vline(xintercept = mean(df_mainland$overlap_weights), color = "#69b3a2") +
  geom_histogram(data = df_barrier, 
                 aes(x = overlap_weights, y = after_stat(count / sum(count))), 
                 fill = "#404080", alpha = 0.7) +
  geom_vline(xintercept = mean(df_barrier$overlap_weights), color = "#404080") +
  theme_bw() 

ggsave("L:/Project-CBRA/output/figures/weights_by_barrier.png", width = 8, height = 5)




################## DIFF-IN-DIFF USING DEVELOPMENT MEASURE #############


df_did <- rbind(est_data %>% select(region_id, 
                                    developed = develpd, 
                                    white = race_white_a,
                                    black = race_black_a,
                                    med_hhincome = med_hhincome_a,
                                    bagrad = ed_ba_gread_a,
                                    overlap_weights),
                est_data %>% select(region_id, 
                                    developed = developed_2018,
                                    white = white_acs,
                                    black = black_acs,
                                    med_hhincome = med_hh_inc_acs,
                                    bagrad = bagrad_acs,
                                    overlap_weights))

df_did$post <- 0
df_did$post[616:1230] <- 1

df_did <- df_did %>% left_join(est_data %>% select(region_id, in_cbrs, barrier),
                               by = "region_id")

reg1 <- felm(developed ~ in_cbrs*post + barrier*post | region_id | 0 | region_id, weight = df_did$overlap_weights, data = df_did) 
reg2 <- felm(white ~ in_cbrs*post + barrier*post | region_id | 0 | region_id, weight = df_did$overlap_weights, data = df_did) 
reg3 <- felm(black ~ in_cbrs*post + barrier*post | region_id | 0 | region_id, weight = df_did$overlap_weights, data = df_did) 
reg4 <- felm(med_hhincome ~ in_cbrs*post + barrier*post | region_id | 0 | region_id, weight = df_did$overlap_weights, data = df_did) 
reg5 <- felm(bagrad ~ in_cbrs*post + barrier*post | region_id | 0 | region_id, weight = df_did$overlap_weights, data = df_did) 

stargazer(reg1, reg2, reg3, reg4, reg5, 
          type = "latex",
          out = "L:/Project-CBRA/output/tables/tab_did.tex",
          keep = c("post", "in_cbrs:post"), omit.stat = "ser")


################## SFHA BUILDING DENSITY REGRESSIONS #############

# regression formula   
form <- as.formula(sprintf(paste0(depvar, " ~ in_cbrs + barrier")))

# weighted least square
wls <- lm_robust(formula = form,
                 data = temp_data,
                 weights = temp_data$overlap_weights)

reg1 <- felm(unit_sfha_density ~ in_cbrs + barrier, weight = est_data$overlap_weights, data = est_data)
reg2 <- felm(spill_sfha_density ~ in_cbrs + barrier, weight = est_data$overlap_weights, data = est_data)
reg3 <- felm(tot_sfha_density ~ in_cbrs + barrier, weight = est_data$overlap_weights, data = est_data)



stargazer(reg1, reg2, reg3, 
          type = "latex",
          dep.var.caption = "SFHA buildings per acre",
          dep.var.labels = c("In unit", "Spillover", "Overall"),
          out = "L:/Project-CBRA/output/tables/tab_sfha_density.tex",
          keep = "in_cbrs", omit.stat = "ser",
          add.lines = list(c("Outcome mean", "0.021", "0.198", "0.181")))
