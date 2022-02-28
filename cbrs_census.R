##################################################
## Project: CBRS
## Author: Sophie Pesek
## Date: Feb 5, 2022
## Script purpose: Process 1980 census data
## Input data: https://data2.nhgis.org/downloads
## Output data:
##################################################

# load/install packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               sf,
               virdis)
sf_use_s2(FALSE)

# set relative working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # sets directory to the current directory
setwd('../') # relative paths to move directory to the root project directory

# read and filter for coastal states
tracts_raw <- st_read("Data/Raw/Census/nhgis0001_shapefile_tl2000_us_tract_1980/US_tract_1980.shp")

census_raw1 <- read_csv("Data/Raw/Census/nhgis0001_csv/nhgis0001_ds104_1980_tract.csv")
census_raw2 <- read_csv("Data/Raw/Census/nhgis0001_csv/nhgis0001_ds107_1980_tract.csv") %>%
  dplyr::select(., -c(2:21))
census_raw3 <- read_csv("Data/Raw/Census/nhgis0001_csv/nhgis0001_ds116_1980_tract.csv") %>%
  dplyr::select(., -c(2:14))

tracts <- tracts_raw %>% 
  filter(NHGISST %in% c("010", "090", "100", "120", "130", "220", "230", "240", "250", 
                        "280", "340", "360", "370", "440", "450", "480", "510"))

joined_tracts <- tracts %>% left_join(census_raw1, by = "GISJOIN") %>%
  left_join(census_raw2, by = "GISJOIN") %>%
  left_join(census_raw3, by = "GISJOIN") %>%
  rowwise() %>% 
  mutate(race_amerind = sum(c_across(C9D003:C9D005), na.rm = T)) %>%
  mutate(race_asian_pi = sum(c_across(C9D006:C9D014), na.rm = T)) %>%
  mutate(hispanic = sum(c_across(C9E002:C9E005), na.rm = T)) %>%
  mutate(labor_force_employed = sum(DHX001, DHX002, DHX005, DHX006, na.rm = T)) %>%
  mutate(labor_force_unemployed = sum(DHX003, DHX007, na.rm = T)) %>%
  mutate(not_in_labor_force = sum(DHX004, DHX008, na.rm = T)) %>%
  dplyr::select(., -c(C9D003:C9D014, C9E002:C9E005, DHX001:DHX004, DHX005:DHX008)) %>%
  rename(race_white = C9D001,
         race_black = C9D002,
         race_other = C9D015,
         not_hispanic = C9E001,
         median_home_value = C8J001,
         ed_elementary = DHM001,
         ed_less_hs = DHM002,
         ed_hs = DHM003,
         ed_less_college = DHM004,
         ed_college = DHM005,
         median_hh_income = DIE001,
         per_cap_income = DIZ001,
         income_above_pov = DI8001,
         income_below_pov = DI8002,
         population = C6W001)

st_write(joined_tracts, "Data/Int/census.shp")

joined_tracts <- joined_tracts %>%
  st_drop_geometry(.)
write_csv(joined_tracts, "Data/Int/census.csv") # write csv
