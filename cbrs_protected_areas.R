##################################################
## Project: CBRS
## Author: Sophie Pesek
## Date: Feb 28, 2022
## Script purpose: Process Protected Areas Database of the United States (PAD-US)
## Input data: https://www.sciencebase.gov/catalog/item/602597f7d34eb12031138e15 
## Output data: Int/protected_areas.csv
##################################################

# load/install packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               sf,
               viridis,
               gdalUtilities)
sf_use_s2(FALSE)
sqmeters2acres = 0.000247105 # constant to convert square meters to acres

# set relative working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # sets directory to the current directory
setwd('../') # relative paths to move directory to the root project directory

# create function to ensure uoutputs are multipolygons (https://gis.stackexchange.com/questions/389814/r-st-centroid-geos-error-unknown-wkb-type-12)
ensure_multipolygons <- function(X) {
  tmp1 <- tempfile(fileext = ".gpkg")
  tmp2 <- tempfile(fileext = ".gpkg")
  st_write(X, tmp1)
  ogr2ogr(tmp1, tmp2, f = "GPKG", nlt = "MULTIPOLYGON")
  Y <- st_read(tmp2)
  st_sf(st_drop_geometry(X), geom = st_geometry(Y))
}

# read cbrs data
cbrs_raw <- st_read("Data/Raw/CBRS-Boundaries/CBRS_Prohibitions_03122019.shp") # read raw data
cbrs <- cbrs_raw %>%
  filter(SU_Date == "1982-10-18") %>% # filter to 1982
  st_transform(., 5070) # transform to albers equal area crs
cbrs <- cbrs %>% mutate(area = st_area(.)) # calculate cbrs parcel areas

# read pad-us data
file_list <- list.files("Data/Raw/PAD-US/states", pattern = "*gdb", full.names = TRUE) # create list of all .gdb files
shapefile_list <- lapply(file_list, st_read) # read all .gdb files
raw_padus <- rbind(shapefile_list[[1]], shapefile_list[[2]], shapefile_list[[3]], shapefile_list[[4]], shapefile_list[[5]],
                   shapefile_list[[6]], shapefile_list[[7]], shapefile_list[[8]], shapefile_list[[9]], shapefile_list[[10]],
                   shapefile_list[[11]], shapefile_list[[12]], shapefile_list[[13]], shapefile_list[[14]], shapefile_list[[15]],
                   shapefile_list[[16]], shapefile_list[[17]], shapefile_list[[18]]) # bind all 18 coastal states together
padus <- st_transform(raw_padus, st_crs(cbrs)) # transform to same coordinate reference system as cbrs
all_padus <- ensure_multipolygons(padus) # ensure multipolygons
padus_1982 <- all_padus %>%
  filter(Date_Est <= 1982 | is.na(Date_Est)) # crop to pad-us areas established before 1982 or with no date
padus_new <- all_padus %>%
  filter(Date_Est > 1982) # crop to pad-us areas established since 1982

# calculate intersection of cbrs and pad-us data
int <- st_intersection(padus_1982, cbrs)
int <- int %>% mutate(protected_area = st_area(.)) # calculate protected area in each cbrs parcel
int <- int %>% group_by(OBJECTID, GAP_Sts, Date_Est, area) %>%
  summarise(total_protected_area = sum(protected_area)) # sum total number of sq meters in each gap status group

# join intersection with original cbrs data
join <- st_join(cbrs, int) # join to get all cbrs parcels back
join$pct_protected <- as.numeric((join$total_protected_area*sqmeters2acres)/(join$area.x*sqmeters2acres)) # calculate % of acres protected under each gap status

# plot data with percent
ggplot() +  
  geom_sf(data = join, lwd = 0.05, aes(fill = pct_protected)) + # thin line weight to better see density fill
  scale_fill_viridis() + theme_bw()
ggsave("Code/Figures/percent_protected_areas.png")

# trim columns and save as csv
out <- join %>%
  mutate(OBJECTID = OBJECTID.x,
         gap_status = GAP_Sts,
         pct_protected = ifelse(is.na(pct_protected), 0, pct_protected)) %>%  # make density 0 if no part is protected
  select(OBJECTID, pct_protected, GAP_Sts) %>%
  st_drop_geometry(.)
write_csv(out, "Data/Int/protected_area1.csv") # write csv

