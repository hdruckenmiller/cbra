##################################################
## Project: CBRS
## Author: Sophie Pesek
## Date: Feb 10, 2022
## Script purpose: Process National Bridge Inventory
## Input data: https://www.fhwa.dot.gov/bridge/nbi.cfm
## Output data: Int/bridge_density.csv
##################################################

# load/install packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               sf,
               viridis)
sf_use_s2(FALSE)
sqmeters2acres = 0.000247105 # constant to convert square meters to acres

# set relative working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # sets directory to the current directory
setwd('../') # relative paths to move directory to the root project directory

# read  cbrs data
cbrs_raw <- st_read("Data/Raw/CBRS-Boundaries/CBRS_Prohibitions_03122019.shp")
cbrs <- cbrs_raw %>% filter(SU_Date == "1982-10-18")  %>% # filter to 1982
  st_transform(., 5070) # transform to albers equal area crs
cbrs <- cbrs %>% mutate(area = st_area(.)) # calculate cbrs parcel areas

# read g roads data
nbi_raw <- read.csv("Data/Raw/Bridge Inventory/fluna_991992-20160919110712.csv")
nbi <- nbi_raw %>%
  filter(YEAR_BUILT_027 <= 1982 & # filter to bridges build before 1982
           LAT_016 != 0 & LONG_017 != 0 & !is.na(LAT_016) & !is.na(LONG_017) & # drop 0 and NA latitude and longitudes
           LAT_016 < 69200000 & LAT_016 > 24000000 &LONG_017 < 126000000 & LONG_017 > 66000000) %>% # trim to the extent of the continental US
  mutate(Latitude = LAT_016/1000000,
         Longitude = -LONG_017/1000000) %>% # transform lat and lon to proper magnitude and direction
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) # set as sf lat and lon with wgs84 projection
nbi <- st_transform(nbi, st_crs(cbrs)) # transform to same coordinate reference system as cbrs

# calculate intersection of cbrs and bridge data
int <- st_intersection(nbi, cbrs) %>%
  count(OBJECTID) # add count of bridges in each cbrs

# join intersection with original cbrs data
out <- st_join(cbrs, int, by = "OBJECTID") # join to get all cbrs parcels back (even those with no bridges)
out$bridge_density <- as.numeric(out$n/(out$area*sqmeters2acres)) # calculate density in bridges per acre

# plot data
ggplot() +  
  geom_sf(data = out, lwd = 0.05, aes(fill = bridge_density)) + # thin line weight to better see density fill
  scale_fill_viridis() + theme_bw()
ggsave("Code/Figures/bridge_density.png")

# trim columns and save as csv
out <- out %>%
  mutate(OBJECTID = OBJECTID.x,
         bridge_count = n,
         bridge_density = ifelse(is.na(bridge_density), 0, bridge_density)) %>%  # make density 0 if no bridges are present in the area
  select(OBJECTID, bridge_count, bridge_density) %>%
  st_drop_geometry(.)
write_csv(out, "Data/Int/bridge_density.csv") # write csv

