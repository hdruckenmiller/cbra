##################################################
## Project: CBRS
## Author: Sophie Pesek
## Date: Feb 28, 2022
## Script purpose: Process Protected Areas Database of the United States (PAD-US)
## Input data: https://www.sciencebase.gov/catalog/item/602597f7d34eb12031138e15 
## Output data: Int/protected_areas_1982.csv (dates <= 1982), protected_areas_new.csv (dates > 1982)
##################################################

# load/install packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               sf,
               viridis,
               gdalUtilities,
               ggplot2)
sf_use_s2(FALSE)
sqmeters2acres = 0.000247105 # constant to convert square meters to acres
meters2miles = 0.000621371 # constant to convert meters to mile

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
protected_all <- rbind(shapefile_list[[1]], shapefile_list[[2]], shapefile_list[[3]], shapefile_list[[4]], shapefile_list[[5]],
                   shapefile_list[[6]], shapefile_list[[7]], shapefile_list[[8]], shapefile_list[[9]], shapefile_list[[10]],
                   shapefile_list[[11]], shapefile_list[[12]], shapefile_list[[13]], shapefile_list[[14]], shapefile_list[[15]],
                   shapefile_list[[16]], shapefile_list[[17]], shapefile_list[[18]]) # bind all 18 coastal states together
padus <- st_transform(raw_padus, st_crs(cbrs)) # transform to same coordinate reference system as cbrs
all_padus <- ensure_multipolygons(padus) # ensure multipolygons
all_padus$Date_Missing <- with(all_padus, ifelse(is.na(as.numeric(Date_Est)), "Yes", "No"))
non_missing <- all_padus %>% filter(Date_Missing == "No")

padus_1982 <- non_missing %>%
  filter(Date_Est <= 1982) # crop to pad-us areas established before 1982 or with no date
# padus_new <- non_missing %>%
#   filter(Date_Est > 1982) # crop to pad-us areas established since 1982
all_padus <- non_missing

# calculate intersection of cbrs and pad-us data for parcels established before 1982
int_1982 <- st_intersection(padus_1982, cbrs)
int_1982 <- int_1982 %>% mutate(protected_area = st_area(.)) # calculate protected area in each cbrs parcel
int_1982 <- int_1982 %>% group_by(OBJECTID, area) %>%
  summarise(total_protected_area = sum(protected_area)) # sum total number of sq meters in each gap status group

# calculate intersection of cbrs and pad-us data for parcels established since 1982
int <- st_intersection(all_padus, cbrs)
int <- int %>% mutate(protected_area = st_area(.)) # calculate protected area in each cbrs parcel
int <- int %>% group_by(OBJECTID, area) %>%
  summarise(total_protected_area = sum(protected_area)) # sum total number of sq meters in each gap status group

# join intersection with original cbrs data for parcels before 1982
join_1982 <- st_join(cbrs, int_1982) # join to get all cbrs parcels back
join_1982$pct_protected <- as.numeric((join_1982$total_protected_area*sqmeters2acres)/(join_1982$area.x*sqmeters2acres)) # calculate % of acres protected under each gap status

# join intersection with original cbrs data for parcels since 1982
join <- st_join(cbrs, int) # join to get all cbrs parcels back
join$pct_protected <- as.numeric((join$total_protected_area*sqmeters2acres)/(join$area.x*sqmeters2acres)) # calculate % of acres protected under each gap status

## plot data with percent
# ggplot() +  
#   geom_sf(data = join_1982, lwd = 0.05, aes(fill = pct_protected)) + # thin line weight to better see density fill
#   scale_fill_viridis() + theme_bw()
# ggsave("Code/Figures/percent_protected_areas.png")

# trim columns and save as csv for parcels from before 1982
out_1982 <- join_1982 %>%
  mutate(OBJECTID = OBJECTID.x,
         pct_protected = ifelse(is.na(pct_protected), 0, pct_protected)) %>%  # make density 0 if no part is protected
  select(OBJECTID, pct_protected) %>%
  st_drop_geometry(.)

# trim columns and save as csv for areas since 1982
out <- join %>%
  mutate(OBJECTID = OBJECTID.x,
         pct_protected = ifelse(is.na(pct_protected), 0, pct_protected)) %>%  # make density 0 if no part is protected
  select(OBJECTID, pct_protected) %>%
  st_drop_geometry(.)

# filter to only cbrs with no protection (1982) -------
unprotec_cbrs <- join_1982 %>%
  mutate(OBJECTID = OBJECTID.x,
         pct_protected = ifelse(is.na(pct_protected), 0, pct_protected)) %>%  # make density 0 if no part is protected
  select(OBJECTID, pct_protected) %>%
  group_by(OBJECTID) %>%
  summarize(any_protec = sum(pct_protected)) %>%
  filter(any_protec == 0)

# find the nearest tract
nearest <- st_nearest_feature(unprotec_cbrs, padus_1982)
distance <- st_distance(unprotec_cbrs, padus_1982[nearest,], by_element = TRUE)

# store the results
cbrs_matched <- unprotec_cbrs %>% 
  cbind(as.data.frame(nearest)) %>% 
  cbind(as.data.frame(distance))

padus_1982_1 <- padus_1982 %>% 
  st_drop_geometry() %>%
  mutate(nearest = row_number()) %>%
  select(nearest)

cbrs_matched1 <- cbrs_matched %>% left_join(padus_1982_1, by = "nearest")
cbrs_matched1 <- cbrs_matched1 %>%
  mutate(distance = as.numeric(distance)*meters2miles) %>%
  st_drop_geometry(.)

out_1982_1 <- left_join(out_1982, cbrs_matched1, by = "OBJECTID")
out_1982_1 <- out_1982_1 %>%
  mutate(distance_to_protected_area = ifelse(is.na(distance), 0, distance)) %>%  # make distance 0 if overlapping with protected area
  select(OBJECTID, pct_protected, distance_to_protected_area)
write_csv(out_1982_1, "Data/Int/protected_area_1982.csv") # write csv

# filter to only cbrs with no protection (all years) -------
unprotec_cbrs <- join %>%
  mutate(OBJECTID = OBJECTID.x,
         pct_protected = ifelse(is.na(pct_protected), 0, pct_protected)) %>%  # make density 0 if no part is protected
  select(OBJECTID, pct_protected) %>%
  group_by(OBJECTID) %>%
  summarize(any_protec = sum(pct_protected)) %>%
  filter(any_protec == 0)

# find the nearest tract
nearest <- st_nearest_feature(unprotec_cbrs, all_padus)
distance <- st_distance(unprotec_cbrs, all_padus[nearest,], by_element = TRUE)

# store the results
cbrs_matched <- unprotec_cbrs %>% 
  cbind(as.data.frame(nearest)) %>% 
  cbind(as.data.frame(distance))

all_padus_1 <- all_padus %>% 
  st_drop_geometry() %>%
  mutate(nearest = row_number()) %>%
  select(nearest)

cbrs_matched1 <- cbrs_matched %>% left_join(all_padus_1, by = "nearest")
cbrs_matched1 <- cbrs_matched1 %>%
  mutate(distance = as.numeric(distance)*meters2miles) %>%
  st_drop_geometry(.)

out_1 <- left_join(out, cbrs_matched1, by = "OBJECTID")
out_1 <- out_1 %>%
  mutate(distance_to_protected_area = ifelse(is.na(distance), 0, distance)) %>%  # make distance 0 if overlapping with protected area
  select(OBJECTID, pct_protected, distance_to_protected_area)
write_csv(out_1, "Data/Int/protected_area_all.csv") # write csv

# # map missing data
# all_padus_new <- all_padus %>% st_drop_geometry() %>%
#   mutate(code = paste(Mang_Type, GAP_Sts, sep = "-"))
# all_padus_new <- all_padus %>% st_drop_geometry() %>%
#   filter(State_Nm == "AL" | State_Nm == "CT" | State_Nm == "DE" | State_Nm == "FL" | State_Nm == "GA" | State_Nm == "LA" | State_Nm == "MA" | State_Nm == "MD" | State_Nm == "ME" |
#            State_Nm == "MS" | State_Nm == "NC" | State_Nm == "NH" | State_Nm == "NJ" | State_Nm == "NY" | State_Nm == "RI" | State_Nm == "SC" | State_Nm == "TX" | State_Nm == "VA")
# 
# all_padus_new <- all_padus %>% complete(Own_Type, Date_Missing, fill = list(count = 0)) %>%
#   group_by(Own_Type, Date_Missing, .drop = FALSE) %>% 
#   summarize(count = n()) %>%  # count records
#   mutate(pct = count[Date_Missing == "No"]/sum(count)) %>%
#   mutate(pct = ifelse(Date_Missing == "Yes", "", pct)) %>%
#   mutate(pct = as.numeric(pct, na.rm = TRUE))
# 
# all_padus_new <- all_padus_new %>%
#   select(Mang_Type, Date_Missing, Shape_Area) %>%
#   group_by(Mang_Type, Date_Missing) %>% 
#   summarize(area = sum(Shape_Area)) %>%  # count records by species
#   mutate(pct = area[Date_Missing == "No"]/sum(area)) %>%
#   mutate(pct = ifelse(Date_Missing == "Yes", "", pct)) %>%
#   mutate(pct = as.numeric(pct, na.rm = TRUE))
# 
# ggplot(all_padus_new, aes(fill=Date_Missing, y=area, x=Mang_Type)) + 
#   geom_bar(position="stack", stat="identity") +
#   scale_fill_manual(values = c("#6FDA87", "#FF5F6B")) +
#   geom_text(aes(label = scales::percent(pct, accuracy = 0.1L)), position = position_stack(vjust = 0.5)) +
#   labs(title="Management of PAD units missing dates by area (coastal states)", x="Manager type", y="Area", fill="Date missing") +
#   theme_classic()
# ggsave("Code/Figures/padus_missing_dates_mang_area.png")
# 
# ggplot(all_padus, aes(fill=Date_Missing, y=frequency(Mang_Type), x=Mang_Type)) + 
#   geom_bar(position="stack", stat="identity") +
#   scale_fill_manual(values = c("#6FDA87", "#FF5F6B")) +
#   labs(title="Distribution of PAD units missing dates", x="Manager type", y="Count", fill="Date missing") +
#   theme_classic()
# ggsave("Code/Figures/padus_missing_dates3.png")
# 
# state <- all_padus %>%
#   filter(Mang_Type == "STAT")
# 
# # map data missing
# ggplot() +  
#   geom_sf(data = state, lwd = 0.05, aes(color = Date_Missing, fill = Date_Missing)) + # thin line weight to better see density fill
#   geom_sf(data = cbrs, lwd = 0.08, fill = NA) +
#   scale_color_manual(values = c("#6FDA87", "#FF5F6B")) +
#   scale_fill_manual(values = c("#6FDA87", "#FF5F6B")) +
#   guides(color = "none") +
#   labs(title="Map of PAD units with missing dates (state)", fill="Date missing") +
#   theme_bw()
# ggsave("Code/Figures/protected_areas_map_state.png", dpi=700)
#
# ggplot(all_padus, aes(x=as.numeric(Date_Est))) + 
#   geom_histogram(binwidth = 5) +
#   xlim(c(1900, 2020)) +
#   #abline(v=1982,col="blue",lwd=3) +
#   geom_vline(aes(xintercept = 1982), color="blue") +
#   labs(title="Histogram of PAD establishment dates", x="Date established", y="Count") +
#   theme_classic()
# ggsave("Code/Figures/padus_date_histogram.png")
