##### Script information

# Paper: Druckenmiller et al (Nature Cliamte Change, 2024) Removing Development Incentives in Risky Areas Promotes Climate Adaptation

# Task: estimate direct effects of CBRS  

# Data inputs
#     data/direct_outcome_data.csv
#     data/synth_weights.csv
#     data/matching_data.csv

# Figure / Table outputs 
#     Table 2

########################################################

# Clear environment, load packages
rm(list=ls())
options(warn=-1)

library(pacman)
pacman::p_load(lfe, microsynth, dplyr, xtable, cobalt, weights)

########################################################


# Load weights 
weights <- read.csv("../data/synth_weights.csv")

# Load outcomes data 
data <- read.csv("../data/direct_outcome_data.csv")

# Load matching data 
df <- read.csv("../data/matching_data.csv")


# Add weights to outcomes data 
data <- merge(data, weights, by = "region_id", all.x = T)
data$synth_weight_ztrax[is.na(data$synth_weight_ztrax)] <- 0

# Outcome variables in full sample 
outcomes <- c("n_buildings_acre",
              "TotalAssessedValue_acre",
              "white_acs", 
              "black_acs", 
              "med_hh_inc_acs", 
              "med_rent_inc_pct_acs",
              "occupied_acs")

# ZTRAX outcomes 
ztrax_outcomes <- c("mean_price", 
                    "LotSizeAcres_MEAN", 
                    "sqfeet_MEAN", 
                    "TotalBedrooms_MEAN")


# Create empty dataframes to fill with estimates 
results <- as.data.frame(matrix(NA, length(outcomes), 7))
colnames(results) <- c("outcome", "estimate", "se", "pvalue", "control_mean", "relative_effect", "N")
results$outcome <- outcomes

ztrax_results <- as.data.frame(matrix(NA, length(ztrax_outcomes), 7))
colnames(ztrax_results) <- c("outcome", "estimate", "se", "pvalue", "control_mean", "relative_effect", "N")
ztrax_results$outcome <- ztrax_outcomes


#### TABLE 2 
#### Estimate effect of CBRS treatment on outcomes using regression with synthetic control weights 


# In the full sample 
for(outcome in outcomes){
  reg <- lm(data = data, 
            formula(paste0(outcome, " ~ in_cbrs + barrier")), weights = synth_weight)
  base <- weighted.mean(data[(data$in_cbrs==0 & !is.na(data$synth_weight)), outcome], data[(data$in_cbrs==0 & !is.na(data$synth_weight)), "synth_weight"], na.rm = T)
  effect <- summary(reg)$coefficients["in_cbrs",1]
  relative <- ((base + effect) - base)/base
  results$estimate[results$outcome==outcome] <- round(effect, 4)
  results$se[results$outcome==outcome] <- round(summary(reg)$coefficients["in_cbrs",2], 4)
  results$pvalue[results$outcome==outcome] <- round(summary(reg)$coefficients["in_cbrs",4], 4)
  results$control_mean[results$outcome==outcome] <- round(base, 4)
  results$relative_effect[results$outcome==outcome] <- round(relative, 4)
  results$N[results$outcome==outcome] <- nobs(reg)
}

# In the ZTRAX sample 
for(outcome in ztrax_outcomes){
  print(outcome)
  reg <- lm(data = data, 
            formula(paste0(outcome, " ~ in_cbrs + barrier")), weights = synth_weight_ztrax)
  base <- weighted.mean(data[data$in_cbrs==0, outcome], data[data$in_cbrs==0, "synth_weight_ztrax"], na.rm = T)
  effect <- summary(reg)$coefficients["in_cbrs",1]
  relative <- ((base + effect) - base)/base
  ztrax_results$estimate[ztrax_results$outcome==outcome] <- round(effect, 4)
  ztrax_results$se[ztrax_results$outcome==outcome] <- round(summary(reg)$coefficients["in_cbrs",2], 4)
  ztrax_results$pvalue[ztrax_results$outcome==outcome] <- round(summary(reg)$coefficients["in_cbrs",4], 4)
  ztrax_results$control_mean[ztrax_results$outcome==outcome] <- round(base, 4)
  ztrax_results$relative_effect[ztrax_results$outcome==outcome] <- round(relative, 4)
  ztrax_results$N[ztrax_results$outcome==outcome] <- nobs(reg)
  
}

# For BUFA (which uses DID instead of cross-sectional estimation)
bufa <- df[(df$year %in% c(1980, 2010)),c("region_id", "year", "bufa_share")]
bufa <- merge(bufa, weights, by = "region_id")
bufa <- merge(bufa, data[,c("region_id", "in_cbrs", "barrier")], by = "region_id")
bufa$treated <- 0
bufa$treated[bufa$in_cbrs==1 & bufa$year==2010] <- 1

reg <- felm(data = bufa,
            formula("bufa_share ~ treated + barrier | in_cbrs + year"), weights = bufa$synth_weight)
base <- weighted.mean(bufa[(bufa$in_cbrs==0 & bufa$year==2010), "bufa_share"], bufa[(bufa$in_cbrs==0 & bufa$year==2010), "synth_weight"], na.rm = T)
effect <- summary(reg)$coefficients["treated",1]
relative <- ((base + effect) - base)/base
bufa_results <- c("bufa_share", 
                  round(effect, 4),
                  round(summary(reg)$coefficients["treated",2], 4),
                  round(summary(reg)$coefficients["treated",4], 4),
                  round(base, 4),
                  round(relative, 4),
                  sum(bufa$synth_weight > 0))

# Combine all direct effect estimates and output Table 2 
out <- rbind(results, ztrax_results, bufa_results)
out
print.xtable(xtable(out), include.rownames = F, file = "../output/Table2.tex")

#################################################
