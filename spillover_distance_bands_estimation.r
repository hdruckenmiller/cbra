# Calculate means in spillover distance bands of CBRS, counterfactuals; estimation for Zillow, NFIP, and census outcomes
# Date last modified: 10/1/2022
# Last modified by: Penny Liao

rm(list = ls())

library(pacman)

p_load(haven, readr, tidyverse, tidymodels, magrittr, janitor, skimr, zoo, ggplot2, RColorBrewer, 
       data.table, gravity, fastDummies, texreg, fixest, sf, ggmap, raster, sp, sf, here, maptools, 
       rgdal, maps, rgeos, lubridate, ggpubr, grid, estimatr)


setwd("L:/Project-CBRA/data")


## Function to run regressions =======================================================================================

est_felm <- function(depvar, label, weights) {
  
  form <- as.formula(sprintf(paste0(depvar, "~ buffer_500 + buffer_1000 + buffer_1500 + buffer_2000 + ", 
                                    "buffer_500:trt + buffer_1000:trt + buffer_1500:trt + buffer_2000:trt + barrier + 0")))
  
  if (weights == "fullweights") {
    model = lm_robust(data = est_data, 
                      formula = form,
                      weights = overlap_weights_full_sample)
  } else {
    model = lm_robust(data = est_data, 
                      formula = form,
                      weights = overlap_weights_zsample)
  }
  summary(model)

  
  fit_mean_avg = as.data.frame(coef(summary(model))[, 1:2])
  
  # no barrier control
  #fit_mean_avg = fit_mean_avg[5:8,]
  
  # with barrier control
  fit_mean_avg = fit_mean_avg[6:9,]
  
  fit_mean_avg = fit_mean_avg %>% mutate(
    Distance = c("500", "1000", "1500", "2000"))
  
  fit_mean_avg = fit_mean_avg %>% rename(
    mean = Estimate,
    sd = `Std. Error`) #%>% mutate(mean = mean/1000,
  #   sd = sd/1000)
  
  
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



## Zillow Main outcomes =======================================================================================

## Read in zillow aggregated data

spillover_zasmt_unit_weighted = read_csv("int/zillow/spillover_zasmt_unit_weighted.csv")
spillover_ztrans_unit_weighted = read_csv("int/zillow/spillover_ztrans_unit_weighted.csv")

# # Merge in weights accounting for barrier identifiers
# fwt_barrier <- read_csv("int/est_weights/weights_fullsample_barrier.csv")
# zwt_barrier <- read_csv("int/est_weights/weights_main_zsample_barrier.csv")
# 
# # region_ids where there is no leading zero in match_sample
# add_leading_zero <- c("154", "1117", "1208", "1219", "1220", "1223", 
#                       "1230", "1231", "1242", "1255", "1267", "1289", "1294", 
#                       "94363", "94367", "95093", "99344")
# 
# fwt_barrier$region_id[fwt_barrier$region_id %in% add_leading_zero] <- paste0("0", fwt_barrier$region_id)
# zwt_barrier$region_id[zwt_barrier$region_id %in% add_leading_zero] <- paste0("0", zwt_barrier$region_id)
# names(fwt_barrier)[2] <- "overlap_weights_full_sample"
# names(zwt_barrier)[2] <- "overlap_weights_zsample"
# 
# spillover_zasmt_unit_weighted <- spillover_zasmt_unit_weighted %>%
#   dplyr::select(-overlap_weights_full_sample, -overlap_weights_zsample) %>%
#   left_join(fwt_barrier, by = c("regin_d" = "region_id")) %>%
#   left_join(zwt_barrier, by = c("regin_d" = "region_id"))


# Merge in barrier island/capes identifiers 
barriers <- read.csv("int/barriers/units_barriers.csv")
barriers$treated <- NULL

spillover_zasmt_unit_weighted <- spillover_zasmt_unit_weighted %>% 
  left_join(barriers, by = c("regin_d" = "region_id"))

spillover_ztrans_unit_weighted <- spillover_ztrans_unit_weighted %>% 
  left_join(barriers, by = c("regin_d" = "region_id"))

rm(barriers)



# no_buildings_per_area (use full sample weights)
est_data <- spillover_zasmt_unit_weighted
no_buildings_per_area <- est_felm("no_buildings_per_area", "a. Buildings per acre", "fullweights")

no_buildings_per_area$data

## avg_sales_price_amount (use zsample weights)
est_data <- spillover_ztrans_unit_weighted %>%
  mutate(avg_sales_price_amount = avg_sales_price_amount/1000)

avg_sales_price_amount <- est_felm("avg_sales_price_amount", "b. Avg. sales price (thou. $)", "zweights")


## total_assessed_value_per_acre (use zsample weights)
est_data <- spillover_zasmt_unit_weighted %>%
  mutate(total_assessed_value_per_acre = total_assessed_value_per_acre/1000)

total_assessed_value_per_acre <- est_felm("total_assessed_value_per_acre", "c. Total assessed value per acre ($1000)", "zweights")



## NFIP Main outcomes =======================================================================================

## Read in NFIP aggregated data
nfip_spillover_0_2000km_aggregate_weighted = read_csv("int/nfip_outcomes/nfip_spillover_0_2000km_aggregate_weighted.csv")


# # Merge in weights accounting for barrier islands
# nfip_spillover_0_2000km_aggregate_weighted <- nfip_spillover_0_2000km_aggregate_weighted %>%
#   left_join(fwt_barrier, by = c("regin_d" = "region_id"))

# Rename weight variable to fit regression function
nfip_spillover_0_2000km_aggregate_weighted = nfip_spillover_0_2000km_aggregate_weighted %>%
  rename(overlap_weights_full_sample = overlap_weights)

# Load inland counterfactual units to drop
drop <- read.csv("int/regionalization/landuse_counterfactuals_ratio3/counterfactual_inland_drop.csv")
names(drop) <- "region_id"
nfip_spillover_0_2000km_aggregate_weighted <- nfip_spillover_0_2000km_aggregate_weighted %>% 
  mutate(region_id = as.numeric(regin_d)) %>%
  filter(!region_id %in% drop$region_id)
rm(drop)

# Merge in barrier island/capes identifiers 
barriers <- read.csv("int/barriers/units_barriers.csv")
barriers$treated <- NULL
nfip_spillover_0_2000km_aggregate_weighted <- nfip_spillover_0_2000km_aggregate_weighted %>%
  left_join(barriers, by = c("regin_d" = "region_id"))
rm(barriers)


## tot_claims_per_acre

est_data <- nfip_spillover_0_2000km_aggregate_weighted

tot_claims_per_acre <- est_felm("tot_claims_per_acre", "d. Tot. claims ($) per acre", "fullweights")


#tot_claims_per_1000_coverage

est_data <- nfip_spillover_0_2000km_aggregate_weighted

tot_claims_per_1000_coverage <- est_felm("tot_claims_per_1000_coverage", "e. NFIP claims per $1000 coverage ($)", "fullweights")



#tot_bldg_per_acre

est_data <- nfip_spillover_0_2000km_aggregate_weighted

tot_bldg_per_acre <- est_felm("tot_bldg_per_acre", "f. Buildings per acre in SFHA", "fullweights")





## Census outcomes ==================================================================================================


## Read in census data

spillover_bands_census2020 = read.csv("int/outcomes/spillover_bands_census2020.csv")


# Extract list of CBRS ID
match_sample <- read.csv("int/match_sample.csv")
match_sample <- match_sample %>% filter(in_cbrs == 1)
cbrs_id <- match_sample$region_id

# Create treatment variable
spillover_bands_census2020 = spillover_bands_census2020 %>% mutate(
  trt = ifelse(regin_d %in% cbrs_id,1,0)
)
rm(match_sample, cbrs_id)


spillover_bands_census2020 =  spillover_bands_census2020 %>% mutate(buffer_500 = ifelse(dstnc_b == 0,1,0),
                                                                    buffer_1000 = ifelse(dstnc_b == 500,1,0),
                                                                    buffer_1500 = ifelse(dstnc_b == 1000,1,0),
                                                                    buffer_2000 = ifelse(dstnc_b == 1500,1,0))


## Combine with full sample weight

## Read in full sample weight updated (plus 0 in front for the first couple of states)

weights_fullsample_updated <- read_csv("int/est_weights/weights_fullsample_updated.csv")

spillover_bands_census2020_weighted = left_join(spillover_bands_census2020, 
                                                weights_fullsample_updated %>% rename(regin_d = region_id),
                                                by = c("regin_d"))

# Rename weight variable to fit regression function
spillover_bands_census2020_weighted = spillover_bands_census2020_weighted %>% 
  rename(overlap_weights_full_sample = overlap_weights)

# # Rename weight variable to fit regression function
# spillover_bands_census2020_weighted = spillover_bands_census2020 %>% 
#   left_join(fwt_barrier, by = c("regin_d" = "region_id"))


# Merge in barrier island/capes identifiers 
barriers <- read.csv("int/barriers/units_barriers.csv")
barriers$treated <- NULL
spillover_bands_census2020_weighted <- spillover_bands_census2020_weighted %>%
  left_join(barriers, by = c("regin_d" = "region_id"))
rm(barriers)


# "occupied_p"

est_data <- spillover_bands_census2020_weighted

occupied_p <- est_felm("occupied_p", "g. Owner-occupied units (%)", "fullweights")


## "white_p"

est_data <- spillover_bands_census2020_weighted

white_p <- est_felm("white_p", "h. White (%)", "fullweights")



## "black_p"

est_data <- spillover_bands_census2020_weighted

black_p <- est_felm("black_p", "h. Black (%)", "fullweights")


## Create PDF ==================================================================================================

pdf("L:/Project-CBRA/output/spatial_lags_model/spatial_lag_result_barrier_ctrl.pdf",         # File name
    width = 8, height = 5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "letter") 


ggarrange( no_buildings_per_area  + rremove("xlab"), avg_sales_price_amount + rremove("xlab"), total_assessed_value_per_acre + rremove("xlab"), 
           tot_claims_per_acre + rremove("xlab"),  tot_claims_per_1000_coverage  + rremove("xlab"), tot_bldg_per_acre + rremove("xlab"), 
           ncol = 3,nrow = 2,
           common.legend = TRUE, legend = "bottom")


dev.off()




## ACS outcomes ==================================================================================================


## Read in ACS data
setwd("L:/Project-CBRA/data")
spillover_bands_acs = read.csv("int/outcomes/spillover_band_acs.csv")

# Extract list of CBRS ID
match_sample <- read.csv("int/match_sample.csv")
match_sample <- match_sample %>% filter(in_cbrs == 1)
cbrs_id <- match_sample$region_id

# Create treatment variable
spillover_bands_acs = spillover_bands_acs %>% mutate(
  trt = ifelse(regin_d %in% cbrs_id,1,0)
)
rm(match_sample, cbrs_id)


spillover_bands_acs =  spillover_bands_acs %>% mutate(buffer_500 = ifelse(dstnc_b == 0,1,0),
                                                      buffer_1000 = ifelse(dstnc_b == 500,1,0),
                                                      buffer_1500 = ifelse(dstnc_b == 1000,1,0),
                                                      buffer_2000 = ifelse(dstnc_b == 1500,1,0))


## Combine with full sample weight

## Read in full sample weight updated (plus 0 in front for the first couple of states)

weights_fullsample_updated <- read_csv("int/est_weights/weights_fullsample_updated.csv")

spillover_bands_acs_weighted = spillover_bands_acs %>%
  left_join(weights_fullsample_updated, by = c("regin_d" = "region_id"))

# Rename weight variable to fit regression function
spillover_bands_acs_weighted = spillover_bands_acs_weighted %>%
  rename(overlap_weights_full_sample = overlap_weights)

# # Merge in weights accounting for barrier islands
# spillover_bands_acs_weighted = spillover_bands_acs %>%
#   left_join(fwt_barrier, by = c("regin_d" = "region_id"))

# Merge in barrier island/capes identifiers 
barriers <- read.csv("int/barriers/units_barriers.csv")
barriers$treated <- NULL
spillover_bands_acs_weighted <- spillover_bands_acs_weighted %>%
  left_join(barriers, by = c("regin_d" = "region_id"))
rm(barriers, weights_fullsample_update, spillover_bands_acs)



## "white"

est_data <- spillover_bands_acs_weighted

white <- est_felm("white", "Share White", "fullweights")


## "black"

est_data <- spillover_bands_acs_weighted

black <- est_felm("black", "Share Black", "fullweights")


## "hispanic"

est_data <- spillover_bands_acs_weighted

hispanic <- est_felm("hispanic", "Share Hispanic", "fullweights")


## "bagrad"

est_data <- spillover_bands_acs_weighted

bagrad <- est_felm("bagrad", "Share college grad", "fullweights")


## "median_hh_income"

est_data <- spillover_bands_acs_weighted

median_hh_income <- est_felm("median_hh_income", "Median household income", "fullweights")



## "own_occupied"

est_data <- spillover_bands_acs_weighted

own_occupied <- est_felm("own_occupied", "Share owner-occupied units", "fullweights")


## "rent_occupied"

est_data <- spillover_bands_acs_weighted

rent_occupied <- est_felm("rent_occupied", "Share renter-occupied units", "fullweights")


## "occupied"

est_data <- spillover_bands_acs_weighted

occupied <- est_felm("occupied", "Share cccupied unit", "fullweights")



## "median_rent"

est_data <- spillover_bands_acs_weighted

median_rent <- est_felm("median_rent", "Median rent", "fullweights")



## "median_rent_inc_pct"
est_data <- spillover_bands_acs_weighted

median_rent_inc_pct <- est_felm("median_rent_inc_pct", "Median rent as % of income", "fullweights")


## Save results in a pdf file

pdf("L:/Project-CBRA/output/spatial_lags_model/spatial_lag_acs_barrier_ctrl.pdf",         # File name
    width = 8, height = 5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "letter") 


ggarrange( white + rremove("xlab"), black + rremove("xlab"), bagrad + rremove("xlab"), 
           median_hh_income  + rremove("xlab"), median_rent + rremove("xlab"), rent_occupied + rremove("xlab"), 
           ncol = 3,nrow = 2,
           common.legend = TRUE, legend = "bottom")


dev.off()


#annotate_figure(figure, bottom = textGrob("Distance to cbra area (m)", rot = 0, vjust = 0.3, gp = gpar(cex = 1.3,fontsize=14)))



## Regression observations


#test = spillover_bands_census2020 %>% group_by(dstnc_b) %>% summarise(count = n())


#Observations: 

#nfip regression: 2519 
#0-500m - 644, 500-1000m - 637, 1000-1500m - 622, 1500-2000m - 615

#zasmt: 2110
#0-500m - 534, 500-1000m - 538, 1000-1500m - 528, 1500-2000m - 510  

#ztrans: 1953
#0-500m - 497, 500-1000m - 501, 1000-1500m - 488, 1500-2000m - 467   

#census: 2466
#0-500m - 628, 500-1000m - 623, 1000-1500m - 612, 1500-2000m - 603


