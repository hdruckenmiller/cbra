##### Script information

# Paper: Druckenmiller et al (Nature Cliamte Change, 2024) Removing Development Incentives in Risky Areas Promotes Climate Adaptation

# Task: generate weights for synthetic control analysis 

# Data inputs
#     data/matching_data.csv
#     data/matching_data_ZTRAX.csv

# Data outputs
#     data/synth_weights.csv
#     output/source_data_fig3a.csv

# Figure / Table outputs 
#     Figure 3a
#     Table 1 
#     Figure ED3
#     Figure ED5 
#     Table ED1

########################################################

# Clear environment, load packages
rm(list=ls())
options(warn=-1)

library(pacman)
pacman::p_load(lfe, microsynth, dplyr, xtable, cobalt, weights, scales)

########################################################


##### Load matching data 
df <- read.csv("../data/matching_data.csv")

##### Pre-process data for synthetic controls 

# Need to divide outcome variables by totals so that scaled correctly for micro_synth package 
df$bufa_share_percent_agg <- (df$bufa_share/length(unique(df$region_id[df$CBRS_unit==1])))
df$spillover_bufa_share_percent_agg <- (df$spillover_bufa_share/length(unique(df$region_id[df$CBRS_unit==1])))

# Normalize variables outside 0-1 range 
df$elevatn_norm <- (df$elevatn - mean(df$elevatn))/sd(df$elevatn)
df$urban_population_mean_norm <- (df$urban_population_mean - mean(df$urban_population_mean))/sd(df$urban_population_mean)
df$pre1982_SUM_norm <- (df$pre1982_SUM - mean(df$pre1982_SUM))/sd(df$pre1982_SUM)
df$med_hhincome_p_norm <- (df$med_hhincome_p - mean(df$med_hhincome_p))/sd(df$med_hhincome_p)

# Define matching variables 
cov.var.norm <- c("wetland", "barren", "treecvr", 
                  "elevatn_norm",  "pre1982_SUM_norm", 
                  "urban_population_mean_norm", "med_hhincome_p_norm",
                  "A_all_share", "V_all_share",
                  "region_NorthAtlantic", "region_SouthAtlantic", "region_GulfCoast")

##### Run synthetic control algorithm  

synth <- microsynth(as.data.frame(df),
                    idvar="region_id", timevar="year", intvar="treated", 
                    match.out=c("bufa_share_percent_agg", "spillover_bufa_share_percent_agg"), 
                    match.covar=cov.var.norm, 
                    end.post = 2010, end.pre = 1980, start.pre = 1960,
                    result.var="bufa_share_percent_agg", 
                    test="lower", maxit = 5000, cal.epsilon = 1e-04,
                    n.cores = min(parallel::detectCores(), 2))

##### Figure 3, panel a  (synthetic controls plots)

pdf("../output/Figure3a.pdf", width = 11, height = 5)
plot_microsynth(synth, start.pre = 1960, end.pre = 1980, end.post = 2010,
                main.tc = "", ylab.tc = "Built-up surface share (%)",
                main.diff = "")
dev.off()

##### Table 1 (balance table)

# Identify variables to include in balance table 
balance.vars.nobufa <- c("wetland", "barren", "treecvr", 
                         "elevatn",  "pre1982_SUM", 
                         "urban_population_mean", "med_hhincome_p",
                         "A_all_share", "V_all_share",
                         "region_NorthAtlantic", "region_SouthAtlantic", "region_GulfCoast")
balance.vars <- c("wetland", "barren", "treecvr", 
                  "elevatn",  "pre1982_SUM", 
                  "urban_population_mean", "med_hhincome_p",
                  "A_all_share", "V_all_share",
                  "region_NorthAtlantic", "region_SouthAtlantic", "region_GulfCoast",
                  "bufa_share_1960", "bufa_share_1970", "bufa_share_1980",
                  "spillover_bufa_share_1960", "spillover_bufa_share_1970", "spillover_bufa_share_1980")

# Function to generate balance table 
balance_table <- function(model, dataframe){
  weights <- as.data.frame(model$w$Weights)
  weights$region_id <- rownames(weights)
  colnames(weights)[1] <- "w"
  balance <- dataframe[dataframe$year==1980, c("region_id", "CBRS_unit",  balance.vars.nobufa)]
  bufa_60 <- dataframe[dataframe$year==1960, c("region_id", "bufa_share", "spillover_bufa_share")]
  colnames(bufa_60) <- c("region_id", "bufa_share_1960", "spillover_bufa_share_1960")
  bufa_70 <- dataframe[dataframe$year==1970, c("region_id", "bufa_share", "spillover_bufa_share")]
  colnames(bufa_70) <- c("region_id", "bufa_share_1970", "spillover_bufa_share_1970")
  bufa_80 <- dataframe[dataframe$year==1980, c("region_id", "bufa_share", "spillover_bufa_share")]
  colnames(bufa_80) <- c("region_id", "bufa_share_1980", "spillover_bufa_share_1980")
  balance <- merge(balance, bufa_60, by = "region_id")
  balance <- merge(balance, bufa_70, by = "region_id")
  balance <- merge(balance, bufa_80, by = "region_id")
  balance <- merge(balance, weights, by = "region_id")
  btable <- as.data.frame(matrix(NA, length(balance.vars), 4))
  colnames(btable) <- c("variable", "cbrs_mean", "control_mean", "SMD")
  row <- 0
  for(var in balance.vars){
    row <- row + 1
    btable$variable[row] <- var
    btable$cbrs_mean[row] <- round(weighted.mean(balance[balance$CBRS_unit==1, var], w = balance[balance$CBRS_unit==1, "w"]), 4)
    btable$control_mean[row] <- round(weighted.mean(balance[balance$CBRS_unit==0, var], w = balance[balance$CBRS_unit==0, "w"]), 4)
  }
  btable$SMD <- round(cobalt::col_w_smd(balance[,balance.vars], balance$CBRS_unit, balance$w), 4)
  btable 
}

# Generate and view balance table 
btable <- balance_table(synth, df)
btable

# Output Table 1 to Latex file 
print.xtable(xtable(btable), include.rownames = F, file = "../output/Table1.tex")



##### ED Figure 5 (placebo test)


# Make base plot (like Figure 3) to add placebo lines to   
control <- as.vector(synth$Plot.Stats$Control["bufa_share_percent_agg",2:7])
treat <- as.vector(synth$Plot.Stats$Treatment["bufa_share_percent_agg",2:7])

diff <- as.vector(synth$Plot.Stats$Treatment["bufa_share_percent_agg",2:7] - synth$Plot.Stats$Control["bufa_share_percent_agg",2:7])
t <- seq(1960,2010,10)

source_data <- as.data.frame(matrix(NA, length(t), 4))
source_data$year <- t
source_data$treat_bufa <- treat
source_data$control_bufa <- control
source_data$difference <- diff
write.csv(source_data, "../output/source_data_fig3a.csv", row.names = F)


pdf("../output/ED_Figure5.pdf", width = 11, height = 5)
par(mfrow=c(1,2))
plot(control ~ t, type = "l", col = "black", lty = 2, lwd = 2,
     ylim = c(0, 0.5), xlab = NA, ylab = "Built-up surface share (%)", yaxt = "na")
axis(2, las = 2)
lines(treat ~ t, col = "red", lwd = 2)


plot(diff ~ t, ylim = c(-0.25, 0.25), type = "l", col = "red", lwd = 2,
     xlab = NA, ylab = "Treatment - Control", yaxt = "n")
axis(2, las = 2)
control_units <- unique(df$region_id[df$CBRS_unit==F])

# Create empty dataframe to fill with placebo test results 
placebo_out <- as.data.frame(matrix(NA, 100, 8))
colnames(placebo_out) <- c("n", seq(1950,2010,10))
placebo_out$n <- 1:100

# Set seed 
set.seed(100)

# Run loop for placebo simulation where 50 control units are randomly assigned to treatment 
for(i in 1:200){
  if (i == 1){count <- 0}
  placebo_treated <- sample(control_units, 50)
  placebo_df <- subset(df, CBRS_unit == F)
  placebo_df$treated[(placebo_df$region_id %in% placebo_treated & placebo_df$year > 1980)] <- T
  psynth <- microsynth(as.data.frame(placebo_df),
                       idvar="region_id", timevar="year", intvar="treated", 
                       match.out=c("bufa_share_percent_agg", "spillover_bufa_share_percent_agg"), 
                       match.covar=cov.var.norm, 
                       end.post = 2010, end.pre = 1980, start.pre = 1960,
                       result.var=c("bufa_share_percent_agg"), 
                       test="lower", maxit = 5000, cal.epsilon = 1e-04,
                       n.cores = min(parallel::detectCores(), 2))
  pdiff <- as.vector(psynth$Plot.Stats$Treatment["bufa_share_percent_agg",] - psynth$Plot.Stats$Control["bufa_share_percent_agg",])
  if(abs(pdiff[3]) < 0.001 & abs(pdiff[4]) < 0.001 & abs(pdiff[2] < 0.001)){
    lines(pdiff[2:7] ~ t[2:7], col = alpha("black", 0.3))
    count <- count + 1
    print(count)
  }
  if(count == 100){break}
}
lines(diff[2:7] ~ t[2:7], col = "red", lwd = 2)
dev.off()


##### Find synthetic control weights for ZTRAX subsample 

# Load ZTRAX subsample (all units with non-missing ZTRAX data)
ztrax_df <- read.csv("../data/matching_data_ZTRAX.csv")

# Pre-process data for microsynth (same as with full sample)
ztrax_df$bufa_share_percent_agg <- (ztrax_df$bufa_share/length(unique(ztrax_df$region_id[ztrax_df$CBRS_unit==1])))
ztrax_df$spillover_bufa_share_percent_agg <- (ztrax_df$spillover_bufa_share/length(unique(ztrax_df$region_id[ztrax_df$CBRS_unit==1])))

ztrax_df$elevatn_norm <- (ztrax_df$elevatn - mean(ztrax_df$elevatn))/sd(ztrax_df$elevatn)
ztrax_df$urban_population_mean_norm <- (ztrax_df$urban_population_mean - mean(ztrax_df$urban_population_mean))/sd(ztrax_df$urban_population_mean)
ztrax_df$pre1982_SUM_norm <- (ztrax_df$pre1982_SUM - mean(ztrax_df$pre1982_SUM))/sd(ztrax_df$pre1982_SUM)
ztrax_df$med_hhincome_p_norm <- (ztrax_df$med_hhincome_p - mean(ztrax_df$med_hhincome_p))/sd(ztrax_df$med_hhincome_p)

# Find synthetic control weights 
ztrax_synth <- microsynth(as.data.frame(ztrax_df),
                          idvar="region_id", timevar="year", intvar="treated", 
                          match.out=c("bufa_share_percent_agg", "spillover_bufa_share_percent_agg"),
                          match.covar=cov.var.norm, 
                          end.post = 2010, end.pre = 1980, start.pre = 1960,
                          result.var="bufa_share_percent_agg", 
                          test="lower", maxit = 10000, cal.epsilon = 1e-04,
                          #jack=TRUE, perm = 100,
                          n.cores = min(parallel::detectCores(), 2))

# Output ED Figure 3b (note that Figure ED 3a is the same as main text Figure 3a)
pdf("../output/ED_Figure3b.pdf", width = 11, height = 5)
plot_microsynth(ztrax_synth, start.pre = 1960, end.pre = 1980, end.post = 2010,
                main.tc = "", ylab.tc = "Built-up surface share (%)",
                main.diff = "")
dev.off()

# Output ED Table 1 (balance in ztrax subsample) 
btable <- balance_table(ztrax_synth, ztrax_df)
btable
print.xtable(xtable(btable), include.rownames = F, file = "../output/Table_ED1.tex")

# Save weights 
weights <- as.data.frame(synth$w$Weights)
weights$region_id <- rownames(weights)
colnames(weights)[1] <- "synth_weight"

ztrax_weights <- as.data.frame(ztrax_synth$w$Weights)
ztrax_weights$region_id <- rownames(ztrax_weights)
colnames(ztrax_weights)[1] <- "synth_weight_ztrax"

w <- merge(weights, ztrax_weights, by = "region_id", all.x = T)
write.csv(w, "../data/synth_weights.csv", row.names = F)

