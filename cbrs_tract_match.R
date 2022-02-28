# Description: this script overlays CBRS units with 1980 census tracts shapefiles to 
#   calculate the distance to nearest tract for each CBRS unit. 
# Created by: Penny Liao
# Date created: 1/26/2022
# Date last modified: 1/26/2022

### SETUP ###

rm(list = ls())

library(tidyverse)
library(sf)
sf_use_s2(FALSE)
library(viridis)
library(patchwork)

setwd("C:\\Users\\yliao\\OneDrive - rff\\Documents\\Projects\\CBRS\\Data\\Raw")


### READ SPATIAL DATA ###

# read CBRS units and filter for original ones
cbrs_raw <- st_read("CBRS-Boundaries/CBRS_Prohibitions_03122019.shp")
cbrs <- cbrs_raw %>% filter(SU_Date == "1982-10-18")

# read and filter for coastal states
tracts_raw <- st_read("Census/nhgis0011_shape/US_tract_1980.shp")
tracts_raw <- st_transform(tracts_raw, st_crs(cbrs))
tracts <- tracts_raw %>% 
  filter(NHGISST %in% c("010", "090", "100", "120", "130", "220", "230", "240", "250", 
                        "280", "340", "360", "370", "440", "450", "480", "510"))

### OVERLAY ###

# find the nearest tract
nearest <- st_nearest_feature(cbrs, tracts)
distance <- st_distance(cbrs, tracts[nearest,], by_element = TRUE)

# store the results
cbrs_matched <- cbrs %>% 
  cbind(as.data.frame(nearest)) %>% 
  cbind(as.data.frame(distance))

tracts1 <- tracts %>% 
  st_drop_geometry() %>%
  mutate(nearest = row_number()) %>%
  select(GISJOIN, nearest)
  
cbrs_matched <- cbrs_matched %>% left_join(tracts1, by = "nearest")
cbrs_matched <- cbrs_matched %>% mutate(distance = as.numeric(distance)/1609.34)

saveRDS(cbrs_matched, "cbrs_nearest_tract.rds")
rm(tracts1)


### SUMMARY ###
cbrs_matched <- readRDS("cbrs_nearest_tract.rds")

summary(cbrs_matched$distance)
cbrs_matched %>% filter(distance == 0) %>% nrow()

# histogram
ggplot(cbrs_matched %>% filter(distance != 0), aes(x = distance)) + 
  geom_histogram(breaks = c(1:50)) 

ggsave("../Meeting notes/tract_dist_hist.png", 
       width = 6, height = 6)

# maps=
centroids <- st_centroid(cbrs_matched)

cbrs_matched <- cbrs_matched %>% mutate(overlap = (distance == 0))
ggplot(centroids) + 
  geom_sf(data = tracts, lwd = 0.1) + 
  geom_sf(aes(color = overlap), size = 0.7, alpha = 0.7)

ggsave("../Meeting notes/tract_overlap.png", 
       width = 8, height = 6)
  
ggplot(centroids) +  
  geom_sf(data = tracts, lwd = 0.1) + 
  geom_sf(aes(color = distance), size = 0.7, alpha = 0.7) +
  scale_color_viridis()

ggsave("../Meeting notes/tract_distance.png", 
       width = 8, height = 6)


### UNIT PLOTS ###

cbrs_matched <- readRDS("cbrs_nearest_tract.rds")

# overlap list
overlap <- st_intersects(cbrs, tracts)

# list of units to plot
units <- c(10, 36, 71, 103)

for (i in 1:length(units)) {
  assign (
    paste0("p", i),
    ggplot(tracts[unlist(overlap[units[i]]),]) +
    geom_sf() +
    geom_sf(data = cbrs[units[i],], fill = "green", alpha = 0.5)
    )
}

(p1 + p2) / (p3 + p4) + plot_annotation(title = "Examples of CBRS units with overlapping tracts")

ggsave("../../Meeting notes/nearest_tract/units_with_overlap.png", 
       width = 8, height = 6)

# non-overlap list
unmatched <- cbrs_matched %>% filter(distance > 0)
unmatched$distance <- round(unmatched$distance, 2)

# create buffer
centroids <- st_centroid(unmatched)
coords <- as.data.frame(st_coordinates(centroids))

# plotting
for (i in 1:nrow(unmatched)) {
  tracts1 <- tracts %>% mutate(nearest = (row_number() == unmatched$nearest[i]))
  
  p1 <- ggplot(unmatched[i,]) + 
    geom_sf(fill = "green", alpha = 0.5) + 
    geom_sf(data = tracts1, aes(color = nearest)) +
    coord_sf(xlim = c(coords[i,1]-2, coords[i,1]+2), ylim = c(coords[i,2]-2, coords[i,2]+2)) +
    ggtitle(paste0("Unit: ", unmatched$Unit[i], ", Distance: ", unmatched$distance[i])) + 
    theme(legend.position="bottom")
  
  p2 <- ggplot(centroids[i,]) +  
    geom_sf(data = tracts, lwd = 0.1) + 
    geom_sf(color = "green") 
  
  p1 + p2
  
  ggsave(paste0("../../Meeting notes/nearest_tract/Unit_", unmatched$Unit[i], ".png"),
         width = 10, height = 5)
}

