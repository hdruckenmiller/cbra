##################################################
## Project: CBRS
## Author: Sophie Pesek
## Date: Feb 26, 2022
## Script purpose: Process National Inventory of Dams
## Input data: https://nid.usace.army.mil/#/
## Output data: Int/dam_density.csv
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

# read cbrs data
cbrs_raw <- st_read("Data/Raw/CBRS-Boundaries/CBRS_Prohibitions_03122019.shp") # read raw data
cbrs <- cbrs_raw %>%
  filter(SU_Date == "1982-10-18") %>% # filter to 1982
  st_transform(., 5070) # transform to albers equal area crs
cbrs <- cbrs %>% mutate(area = st_area(.)) # calculate cbrs parcel areas

# read g roads data
dams_raw <- st_read("Data/Raw/NID/nation.gpkg") # read raw data
dams <- st_transform(dams_raw, st_crs(cbrs)) # transform to same coordinate reference system as cbrs
dams_desc <- read.csv("Data/Raw/NID/nation.csv", skip = 1) # read dam attributes
dams <- left_join(dams, dams_desc, by=c("federalId" = "Federal.ID")) %>% # join dam geolocations with attributes
  filter(Year.Completed <= 1982 | is.na(Year.Completed)) # filter to dams built before 1982 or at an unknown date

# calculate intersection of cbrs and roads data
int <- st_intersection(dams, cbrs)



###
# join intersection with original cbrs data
join <- st_join(cbrs, int, left = TRUE) # join to get all cbrs parcels back (even those with no roads)
out <- group_by(join, OBJECTID.x, area.x) %>%
  summarize(total_length = sum(road_length)) # sum all road lengths within each cbrs area
out$road_density <- as.numeric((out$total_length*meters2miles)/(out$area.x*sqmeters2acres)) # calculate density in miles per acre

ggplot() +  
  geom_sf(data = dams, alpha = 0.2, size = .01, color = "red") + # thin line weight to better see density fill
  geom_sf(data = cbrs, lwd = 0.05) + # thin line weight to better see density fill
  theme_bw()
ggsave("Code/Figures/dam_density.png")


# plot data
ggplot() +  
  geom_sf(data = out, lwd = 0.05, aes(fill = road_density)) + # thin line weight to better see density fill
  scale_fill_viridis() + theme_bw()
ggsave("Code/Figures/dam_density.png")

# trim columns and save as csv
out <- out %>%
  mutate(OBJECTID = OBJECTID.x,
         road_density = ifelse(is.na(road_density), 0, road_density)) %>%  # make density 0 if no roads are present in the area
  select(OBJECTID, road_density) %>%
  st_drop_geometry(.)
write_csv(out, "Data/Int/road_density.csv") # write csv

