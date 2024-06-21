##### Script information

# Paper: Druckenmiller et al (Nature Cliamte Change, 2024) Removing Development Incentives in Risky Areas Promotes Climate Adaptation

# Task: estimate spillover effects 

# Data inputs
#     data/matching_data.csv
#     data/spillover_area_data.csv
#     data/synth_weights.csv


# Data outputs 
#     output/source_data_fig4.csv

# Figure / Table outputs 
#     Figure4.pdf

########################################################

# Clear environment, load packages
rm(list=ls())
options(warn=-1)

library(pacman)
pacman::p_load(estimatr, readr, tidyverse, data.table, ggplot2, egg, ggpubr)

########################################################

#### Function to run spatial lag model regressions and return figure 

est_source <- function(depvar, label, weights) {
  
  form <- as.formula(sprintf(paste0(depvar, "~ buffer_500 + buffer_1000 + buffer_1500 + buffer_2000 + ", 
                                    "buffer_500:trt + buffer_1000:trt + buffer_1500:trt + buffer_2000:trt + barrier + 0 + 
                                    sfha_land_share + padus_land_share + wetland_share + barren_share")))
  
  if (weights == "fullweights") {
    model = lm_robust(data = est_data, 
                      formula = form,
                      weights = synth_weight)
  } else {
    model = lm_robust(data = est_data, 
                      formula = form,
                      weights = synth_weight_ztrax)
  }
  summary(model)
  
  
  fit_mean_avg = as.data.frame(coef(summary(model))[, 1:2])
  
  fit_mean_avg = fit_mean_avg[10:13,]
  
  fit_mean_avg = fit_mean_avg %>% mutate(
    Distance = c("500", "1000", "1500", "2000"))
  
  fit_mean_avg = fit_mean_avg %>% rename(
    mean = Estimate,
    sd = `Std. Error`) 
  
  
  fit_mean_avg$Distance <- factor(fit_mean_avg$Distance, levels=c("500", "1000", "1500", "2000"))
  
  p <- ggplot(fit_mean_avg, aes(x= Distance, y=mean)) +
    geom_errorbar(aes(ymin=mean-sd*1.96, ymax=mean+sd*1.96), width=0, height =0, size = 1.5, color="coral3", alpha = 0.3) +
    geom_errorbar(aes(ymin=mean-sd*1.64, ymax=mean+sd*1.64), width=0, height =0, size = 1.5, color="coral3") +
    geom_line() +
    geom_point(color = "coral3", size = 2) +
    theme_classic() +
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10),
          plot.title = element_text(size = 10,
                                    color = "black",
                                    face = "bold",
                                    vjust = 0.5,
                                    hjust = 0))  +
    ggtitle(label) +
    geom_hline(yintercept = 0, size=1, linetype="dashed") +
    xlab("Distance to cbra area (m)") + 
    ylab("")

  return(p)
  
}

est_data <- read.csv("../data/spillover_area_data.csv")
weights <- read.csv("../data/synth_weights.csv")
colnames(weights)[which(colnames(weights)== "region_id")] <- "regin_d"
est_data <- merge(est_data, weights, by= "regin_d", all.x = T)

### Panel A : development and property values 

# no_buildings_per_area (use full sample weights)
no_buildings_per_area <- est_source("nbuildings_acres", "a. Buildings per acre", "fullweights")

## avg_sales_price_amount (use zsample weights)
est_data <- est_data %>% 
  mutate(avg_sales_price_amount = avg_sales_price_amount/1000)

avg_sales_price_amount <- est_source("avg_sales_price_amount", "b. Avg. sales price ($1000)", "zweights")

## total_assessed_value_per_acre (use zsample weights)
est_data <- est_data %>%
  mutate(total_assessed_value_per_acre = total_assessed_value_per_acre/1000)

total_assessed_value_per_acre <- est_source("total_assessed_value_per_acre", "c. Total assessed value per acre ($1000)", "fullweights")


### Panel B : NFIP claims  

## tot_claims_per_acre
tot_claims_per_acre <- est_source("tot_claims_per_acre", "d. Tot. claims ($) per acre", "fullweights")

#tot_claims_per_1000_coverage
tot_claims_per_1000_coverage <- est_source("tot_claims_per_1000_coverage", "e. NFIP claims per $1000 coverage ($)", "fullweights")

#tot_bldg_per_acre
tot_bldg_per_acre <- est_source("tot_bldg_per_acre", "f. Buildings per acre in SFHA", "fullweights")


### Panel C : Sociodemographics 

## "median_hh_income"
median_hh_income <- est_source("median_hh_income", "Median household income", "fullweights")


## "occupied"
occupied <- est_source("occupied", "Share cccupied unit", "fullweights")

## "white"
white <- est_source("white", "Share White", "fullweights")

## Save results in a pdf file

pdf("../output/Figure4.pdf",         
    width = 8, height = 5, 
    bg = "white",          
    colormodel = "cmyk",    
    paper = "letter") 


ggarrange( no_buildings_per_area + rremove("xlab"), avg_sales_price_amount + rremove("xlab"), total_assessed_value_per_acre + rremove("xlab"), 
           tot_claims_per_acre  + rremove("xlab"), tot_claims_per_1000_coverage + rremove("xlab"), tot_bldg_per_acre + rremove("xlab"), 
           median_hh_income  + rremove("xlab"), occupied + rremove("xlab"), white + rremove("xlab"), 
           ncol = 3,nrow = 3,
           common.legend = TRUE, legend = "bottom")


dev.off()

###########################

# Source data 

#### Function to run spatial lag model regressions and return source data  

est_source <- function(depvar, label, weights) {
  
  form <- as.formula(sprintf(paste0(depvar, "~ buffer_500 + buffer_1000 + buffer_1500 + buffer_2000 + ", 
                                    "buffer_500:trt + buffer_1000:trt + buffer_1500:trt + buffer_2000:trt + barrier + 0 + 
                                    sfha_land_share + padus_land_share + wetland_share + barren_share")))
  
  if (weights == "fullweights") {
    model = lm_robust(data = est_data, 
                      formula = form,
                      weights = synth_weight)
  } else {
    model = lm_robust(data = est_data, 
                      formula = form,
                      weights = synth_weight_ztrax)
  }
  summary(model)
  
  
  fit_mean_avg = as.data.frame(coef(summary(model))[, 1:2])
  
  fit_mean_avg = fit_mean_avg[10:13,]
  
  fit_mean_avg = fit_mean_avg %>% mutate(
    Distance = c("500", "1000", "1500", "2000"))
  
  fit_mean_avg = fit_mean_avg %>% rename(
    mean = Estimate,
    sd = `Std. Error`) 
  
  
  fit_mean_avg$Distance <- factor(fit_mean_avg$Distance, levels=c("500", "1000", "1500", "2000"))
  
  fit_mean_avg$var <- depvar 
  
  control <- subset(est_data, trt==0)
  agg <- aggregate(control[,depvar], by = list(control$dstnc_b), FUN = "mean", na.rm = T)
  colnames(agg) <- c("Distance", "control_mean") 
  agg$Distance <- agg$Distance + 500
  fit_mean_avg <- merge(fit_mean_avg, agg, by = "Distance")
  fit_mean_avg <- arrange(fit_mean_avg, Distance)
  
  return(fit_mean_avg)
  
}

est_data <- read.csv("../data/spillover_area_data.csv")
weights <- read.csv("../data/synth_weights.csv")
colnames(weights)[which(colnames(weights)== "region_id")] <- "regin_d"
est_data <- merge(est_data, weights, by= "regin_d", all.x = T)

### Panel A : development and property values 

# no_buildings_per_area (use full sample weights)
no_buildings_per_area <- est_source("nbuildings_acres", "a. Buildings per acre", "fullweights")

## avg_sales_price_amount (use zsample weights)
est_data <- est_data %>% 
  mutate(avg_sales_price_amount = avg_sales_price_amount/1000)

avg_sales_price_amount <- est_source("avg_sales_price_amount", "b. Avg. sales price ($1000)", "zweights")

## total_assessed_value_per_acre (use zsample weights)
est_data <- est_data %>%
  mutate(total_assessed_value_per_acre = total_assessed_value_per_acre/1000)

total_assessed_value_per_acre <- est_source("total_assessed_value_per_acre", "c. Total assessed value per acre ($1000)", "fullweights")


### Panel B : NFIP claims  

## tot_claims_per_acre
tot_claims_per_acre <- est_source("tot_claims_per_acre", "d. Tot. claims ($) per acre", "fullweights")

#tot_claims_per_1000_coverage
tot_claims_per_1000_coverage <- est_source("tot_claims_per_1000_coverage", "e. NFIP claims per $1000 coverage ($)", "fullweights")

#tot_bldg_per_acre
tot_bldg_per_acre <- est_source("tot_bldg_per_acre", "f. Buildings per acre in SFHA", "fullweights")


### Panel C : Sociodemographics 

## "median_hh_income"
median_hh_income <- est_source("median_hh_income", "Median household income", "fullweights")


## "occupied"
occupied <- est_source("occupied", "Share cccupied unit", "fullweights")

## "white"
white <- est_source("white", "Share White", "fullweights")

## Source data for figure 4
source_data <- rbind(no_buildings_per_area,
                     avg_sales_price_amount,
                     total_assessed_value_per_acre,
                     tot_claims_per_acre,
                     tot_claims_per_1000_coverage,
                     tot_bldg_per_acre,
                     median_hh_income,
                     occupied,
                     white)
source_data <- source_data[,c("var", "Distance", "mean", "sd", "control_mean")]
colnames(source_data) <- c("outcome_variable", "distance_band", "estimate", "se", "control_mean")
source_data$relative_effect <- source_data$estimate/source_data$control_mean
write.csv(source_data, file = "../output/source_data_fig4.csv", row.names = F)
