##### Script information

# Paper: Druckenmiller et al (Nature Cliamte Change, 2024) Removing Development Incentives in Risky Areas Promotes Climate Adaptation

# Task: estimate individual treatment effects, heterogeneity analysis 

# Data inputs
#     data/matching_data.csv
#     data/heterogeneity_analysis_data.csv

# Data outputs
#     output/source_data_fig3b.csv

# Figure / Table outputs 
#     Figure3b.pdf
#     Table_ED2.tex
#     Table_ED3.tex

########################################################

# Clear environment, load packages
rm(list=ls())
options(warn=-1)

library(pacman)
pacman::p_load(lfe, microsynth, dplyr, xtable, cobalt, weights)

########################################################

### Find individual treatment effects 

# Load matching data 
df <- read.csv("../data/matching_data.csv")

# Pre-process data for synthetic controls 

cov.var <- c("wetland", "barren", "treecvr", 
             "elevatn",  
             "urban_population_mean", 
             "region_NorthAtlantic", "region_SouthAtlantic", "region_GulfCoast")

df <- df[complete.cases(df[,c("region_id", "year", "bufa_share", 
                              "treated", "CBRS_unit", cov.var)]),]

df$elevatn_norm <- (df$elevatn - mean(df$elevatn))/sd(df$elevatn)
df$urban_population_mean_norm <- (df$urban_population_mean - mean(df$urban_population_mean))/sd(df$urban_population_mean)

cov.var.norm <- c("wetland", "barren", "treecvr", 
                  "elevatn_norm",  
                  "urban_population_mean_norm", 
                  "region_NorthAtlantic", "region_SouthAtlantic", "region_GulfCoast")

# Empty dataframes to fill with individual treatment effects
treated_units <- unique(df$region_id[df$CBRS_unit == 1])
years <- seq(1950, 2010, 10)

treatment <- as.data.frame(matrix(NA, length(treated_units)*length(years), 4))
colnames(treatment) <- c("region_id", "year", "treatment", "bufa_share")
treatment$region_id <- rep(treated_units, each = length(years))
treatment$year <- rep(years, length(treated_units))
treatment$treatment <- 1

control <- as.data.frame(matrix(NA, length(treated_units)*length(years), 4))
colnames(control) <- c("region_id", "year", "treatment", "bufa_share")
control$region_id <- rep(treated_units, each = length(years))
control$year <- rep(years, length(treated_units))
control$treatment <- 0

# Loop over each CBRS treatment unit to estimate individual treatment effects 
set.seed(100)
for(unit in treated_units){
  print(paste("Finding synthetic control for unit", which(treated_units == unit), "of", length(treated_units)))
  print(unit)
  data <- rbind(subset(df, region_id == unit), subset(df, CBRS_unit==0))
  synth <- microsynth(as.data.frame(data),
                      idvar="region_id", timevar="year", intvar="treated", 
                      match.out="bufa_share", match.covar=cov.var.norm, 
                      end.post = 2010, end.pre = 1980, start.pre = 1970,
                      result.var="bufa_share", 
                      test="lower", maxit = 5000, cal.epsilon = 1e-04,
                      n.cores = min(parallel::detectCores(), 2))
  treatment$bufa_share[(treatment$region_id==unit & treatment$year==1950)] <- synth$Plot.Stats$Treatment[1,1]
  treatment$bufa_share[(treatment$region_id==unit & treatment$year==1960)] <- synth$Plot.Stats$Treatment[1,2]
  treatment$bufa_share[(treatment$region_id==unit & treatment$year==1970)] <- synth$Plot.Stats$Treatment[1,3]
  treatment$bufa_share[(treatment$region_id==unit & treatment$year==1980)] <- synth$Plot.Stats$Treatment[1,4]
  treatment$bufa_share[(treatment$region_id==unit & treatment$year==1990)] <- synth$Plot.Stats$Treatment[1,5]
  treatment$bufa_share[(treatment$region_id==unit & treatment$year==2000)] <- synth$Plot.Stats$Treatment[1,6]
  treatment$bufa_share[(treatment$region_id==unit & treatment$year==2010)] <- synth$Plot.Stats$Treatment[1,7]
  control$bufa_share[(control$region_id==unit & control$year==1950)] <- synth$Plot.Stats$Control[1,1]
  control$bufa_share[(control$region_id==unit & control$year==1960)] <- synth$Plot.Stats$Control[1,2]
  control$bufa_share[(control$region_id==unit & control$year==1970)] <- synth$Plot.Stats$Control[1,3]
  control$bufa_share[(control$region_id==unit & control$year==1980)] <- synth$Plot.Stats$Control[1,4]
  control$bufa_share[(control$region_id==unit & control$year==1990)] <- synth$Plot.Stats$Control[1,5]
  control$bufa_share[(control$region_id==unit & control$year==2000)] <- synth$Plot.Stats$Control[1,6]
  control$bufa_share[(control$region_id==unit & control$year==2010)] <- synth$Plot.Stats$Control[1,7]
}

### Find units for which we were able to identify an individual synthetic control 
### Meaning units where we were able to reasonable match development pre-trends 

# Calculate percent difference between built-up surface areas in pre-policy periods (1980, 1970)
pchange1980 <- (treatment$bufa_share[treatment$year==1980] - control$bufa_share[control$year==1980])/control$bufa_share[control$year==1980]
pchange1980[is.na(pchange1980)] <- 0
pchange1980[pchange1980 > 1] <- 1

pchange1970 <- (treatment$bufa_share[treatment$year==1970] - control$bufa_share[control$year==1970])/control$bufa_share[control$year==1970]
pchange1970[is.na(pchange1970)] <- 0
pchange1970[pchange1970 > 1] <- 1

# Count as converted if percent difference is less than 10 in both pre-periods 
converged <- (abs(pchange1980) < 0.1) & (abs(pchange1970) < 0.1)

# How many converged? 
print(paste("Indvidual synthetic controls found for",  sum(converged), "units"))

# Subset to units with individual synthetic controls that converged 
treatment <- treatment[treatment$region_id %in% treated_units[converged],]
control <- control[control$region_id %in% treated_units[converged],]

# Drop units with missing built up surface data 
treatment <- treatment[!is.na(treatment$bufa_share),]
control <- control[!is.na(control$bufa_share),]

# Calculate individual treatment effects (absolute)
diff <- treatment$bufa_share[treatment$year==2010] - control$bufa_share[control$year==2010]

# Calculate relative effects 
pchange <- (treatment$bufa_share[treatment$year==2010] - control$bufa_share[control$year==2010])/control$bufa_share[control$year==2010]

# Top code relative effects at 100%
pchange_capped <- pchange
pchange_capped[pchange_capped > 1] <- 1 
#diff_Save <- diff 

# Report mean absolute and relative effects 
mean_effect_bufa <- mean(diff)
mean_bufa_pchange <- (mean(treatment$bufa_share[treatment$year==2010]) - mean(control$bufa_share[control$year==2010]))/mean(control$bufa_share[control$year==2010])
print(paste("Mean treatment effect for bufa share", round(mean_effect_bufa*100, 3),
            ", Mean percent change in bufa share", round(mean_bufa_pchange, 3)))

#### Figure 3b, individual treatment effect plots 
pdf("../output/Figure3b.pdf", width = 11, height = 7)
par(mfrow=c(1,2))
hist(diff, xlab = "Effect on built-up surface (p.p.)", breaks = seq(-1,1, .1), main = NA,
     yaxt = "n")
axis(2, las = 2)
abline(v = 0, lty = 2, col = "red")
hist(pchange_capped*100, xlab = "Relative effect (percent change)", 
     breaks = seq(-100,100, 10), main = NA,
     yaxt = "n")
axis(2, las =2)
abline(v = 0, lty = 2, col = "red")
abline(v = quantile(pchange_capped, 0.25)*100, lty = 2, col = "red")
dev.off()

# Save individual treatment effects 
ind_matches <- unique(treatment$region_id)
het <- treatment %>% filter(year == 1980)
het$pchange <- pchange_capped
het$estimate <- (treatment$bufa_share[treatment$year==2010] - control$bufa_share[control$year==2010])
het <- het[,c("region_id", "bufa_share", "estimate", "pchange")]
colnames(het) <- c("region_id", "bufa_share_1980", "estimate", "pchange")
write.csv(het, "../output/source_data_fig3b.csv", row.names = F)


########

data <- read.csv("../data/heterogeneity_analysis_data.csv")

data <- merge(het, data, by = "region_id")

# Variables to compare in heterogeneity analysis 
vars <- c("region_id", "bufa_share_1980", "wetland", "barren", "elevatn",
          "dst_t_c", "state", "region_NorthAtlantic", "region_SouthAtlantic",
          "region_GulfCoast", "med_hhincome_p", "land_acres", "pop_density",
          "urban_population_mean", "barrier", "mean_price")

# Classify units as effective vs. noneffective 
data$effective <- 0
data$effective[data$pchange < quantile(data$pchange, 0.25)] <- 1

data$non_effective <- 0
data$non_effective[data$pchange > 0] <- 1

### Table ED3 

# Compare effectiveness by state 
data$count <- 1
state <- aggregate(data[,c("estimate", "pchange")], by = list(data$state), FUN = "mean")
colnames(state) <- c("state", "estimate", "pchange")
state_count <- aggregate(data$count, by = list(data$state), FUN = "sum")
colnames(state_count) <- c("state", "n_units")

state <- merge(state, state_count, by = "state")
state
print.xtable(xtable(state), include.rownames = F, file = "../output/Table_ED3.tex")


### Table ED2 

# Compare characteristics of effective vs. noneffective units 
variables <- c("bufa_share_1980", "wetland", "barren", 
               "elevatn", "dst_t_c", "barrier", "log_acres",
               "med_hhincome_p", "mean_price", "log_pop_density",
               "urban_population_mean", "LCV_1982", "reliance")
data$log_acres <- log(data$land_acres + 1)
data$log_pop_density <- log(data$pop_density + 1)

# Subset dataframe to only effective and noneffective units 
clan <- subset(data, effective == 1 | non_effective == 1)

# Empty dataframe to fill with results 
compare <- as.data.frame(matrix(NA, length(variables), 4))
colnames(compare) <- c("variable", "effective_mean", "noneffective_mean", "pvalue")
compare$variable <- variables 

# Loop over variables to fill dataframe with comparisons 
for(var in variables){
  print(var)
  t <- t.test(formula(paste0(var, "~ effective")), data = clan)
  compare$effective_mean[compare$variable==var] <- t$estimate[2]
  compare$noneffective_mean[compare$variable==var] <- t$estimate[1]
  compare$pvalue[compare$variable==var] <- t$p.value
}
compare$difference <- round(compare$effective_mean - compare$noneffective_mean, 5)
compare$pvalue <- round(compare$pvalue, 3)
compare
print.xtable(xtable(compare), include.rownames = F, file = "../output/Table_ED2.tex")

