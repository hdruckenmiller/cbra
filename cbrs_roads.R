##################################################
## Project: CBRS
## Author: Sophie Pesek
## Date: Feb 10, 2022
## Script purpose: Process Global Roads Open Access Data Set
## Input data: https://sedac.ciesin.columbia.edu/data/set/groads-global-roads-open-access-v1
## Output data: Int/road_density.csv
##################################################

# load/install packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               sf,
               viridis)
sf_use_s2(FALSE)
sqmeters2acres = 0.000247105 # constant to convert square meters to acres
meters2miles = 0.000621371 # constant to convert meters to mile

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
roads_raw <- st_read("Data/Raw/Global Roads Open Access Data Set/groads-v1-americas-shp/gROADS-v1-americas.shp") # read raw data
roads <- st_transform(roads_raw, st_crs(cbrs)) # transform to same coordinate reference system as cbrs

# calculate intersection of cbrs and roads data
int <- st_intersection(roads, cbrs)
int <- int %>% mutate(road_length = st_length(.)) # calculate length of roads in each intersection in meters

# join intersection with original cbrs data
join <- st_join(cbrs, int, left = TRUE) # join to get all cbrs parcels back (even those with no roads)
out <- group_by(join, OBJECTID.x, area.x) %>%
  summarize(total_length = sum(road_length)) # sum all road lengths within each cbrs area
out$road_density <- as.numeric((out$total_length*meters2miles)/(out$area.x*sqmeters2acres)) # calculate density in miles per acre

# plot data
ggplot() +  
  geom_sf(data = out, lwd = 0.05, aes(fill = road_density)) + # thin line weight to better see density fill
  scale_fill_viridis() + theme_bw()
ggsave("Code/Figures/road_density.png")

# trim columns and save as csv
out <- out %>%
  mutate(OBJECTID = OBJECTID.x,
         road_density = ifelse(is.na(road_density), 0, road_density)) %>%  # make density 0 if no roads are present in the area
  select(OBJECTID, road_density) %>%
  st_drop_geometry(.)
write_csv(out, "Data/Int/road_density.csv") # write csv
















#function
get_density <- function(x) {
  
  require(tidyverse)
  require(lubridate)
  require(sf)
  
  grids <- cbrs
  x <- 3
  lines <- roads
  
  sub_grids <- grids %>%
    dplyr::filter(OBJECTID == x)
  
  single_lines_objectid <- lines %>%
    #dplyr::filter(OBJECTID == x) %>%
    sf::st_intersection(., sub_grids) %>%
    dplyr::select(OBJECTID, Unit, CBRS_Type, FI_Date, SU_Date, QC, SHAPE_Leng, SHAPE_Area, ROADID, Shape_Leng) %>%
    dplyr::mutate(length_line = st_length(.),
                  length_line = ifelse(is.na(length_line), 0, length_line))
  
  sub_grids <- sub_grids %>%
    sf::st_join(., single_lines_objectid, join = st_intersects) %>%
    dplyr::mutate(OBJECTID = OBJECTID.x, area = SHAPE_Area.y) %>%
    dplyr::group_by(OBJECTID, area) %>%
    dplyr::summarize(length_line = sum(length_line)) %>%
    dplyr::mutate(density = length_line/area) %>%
    # dplyr::mutate(pixel_area = as.numeric(st_area(geom)),
    #               density = length_line/pixel_area) %>%
    dplyr::select(OBJECTID, length_line, area, density)
  return(sub_grids)
}


sapply(1:3, function(j) total.table(scenarios[j]))

new_cbrs <- st_join(cbrs, sub_grids) %>%
  mutate(OBJECTID = OBJECTID.x) %>% select(., -c(OBJECTID.x, OBJECTID.y))




roads_int <- roads %>% st_join(., cbrs, join = st_intersects) %>%
  select(ROADID, LENGTH_KM, Shape_Leng)

cbrs_list <- list(cbrs$OBJECTID)   
  
my_list <- list(df) 

df <- cbrs$OBJECTID
my_list <- list(df) 


for(i in 1:nrow(cbrs)) {             # Using for-loop to add columns to list
  df[i, ] <- cbrs[i, 1]
}

  cbrs %>%
  split(., .$OBJECTID)

sfInit(parallel = TRUE, cpus = num_cores)
sfExport('transmission_lines_hex')
sfSource('src/functions/helper_functions.R')

road_density <- lapply(my_list,
                       function (input_list) {
                         require(tidyverse)
                         require(magrittr)
                         require(lubridate)
                         require(lubridate)
                         require(sf)
                         input_list = my_list
                                       
                         sub_grid <- bind_cols(input_list)
                         unique_ids <- unique(sub_grid$...1)
                                       
                         got_density <- lapply(my_list,
                                               FUN = get_density)
                                       
                         return(got_density)
                         }
)







int = st_intersection(roads, cbrs)


#cropped_roads <- gIntersection(cbrs, roads)

cbrs$area <- st_area(cbrs$geometry) # calculate shape area (is this more or less accurate than the shape area column?)

density_calc <- cbrs %>% st_intersection(roads, cbrs) %>% group_by(area) %>% summarise(road_density = sum(SHAPE_Leng)/sum(area))

road_intersect <- st_intersection(roads, cbrs) %>% # Find the intersections, which should all be points or multilines
  mutate(len_m = st_length(geom)) #%>% # Find the length of each line
  group_by(geometry) %>% # Here you need to insert all the columns from your shapes
  summarize(len_m = sum(len_m))
  
ggplot() + 
  #geom_sf(data = out, lwd = 0.1) + 
  geom_sf(data = cbrs)
  
ggplot(roads) + 
  geom_sf(data = roads, lwd = 0.1) + 
  geom_sf(data = cbrs, color = "blue", alpha = 0.7)

ggsave()
  #geom_sf(aes(color = overlap), size = 0.7, alpha = 0.7)
