library(sf)
library(tidyverse)
#library(magrittr)
#library(ggpubr)
#library(ggplot2)

here::i_am('Fishing_Weights_MCMC.Rproj')

##############################################################
# This script is used to process the raw geospatial data     #
# related to identifying popular fishing lakes in California.#
# These data frames are used in several analyses throughout  #
# the profile and this script exports the data frames for    #
# these analyses.                                            #
##############################################################

## Load boundary data

# load raw UCB CA state boundary
CA.bound = 
  read_sf(here::here('Data/Geospatial/raw/tl_2023_us_state/', 
                     'tl_2023_us_state.shp')) %>%
  filter(NAME == 'California') %>%
  # Reproject to NAD_1983_California_Teale_Albers (meters)
  st_transform(., crs = 3310)

# export California boundary
save(CA.bound, 
     file = here::here('Data/Geospatial/', 'tl_2023_California.Rda'))


## Load and process fishing locations data

# load raw CDFW fishing spots
Fishing.Locs = 
  read_sf(here::here('Data/Geospatial/raw', 
                     'CAfishinglocations.shp')) %>%
  # Reproject to NAD_1983_California_Teale_Albers (meters)
  st_transform(., crs = 3310) %>%
  st_zm(.) %>% # remove z dimension
  # select only locations that have been stocked during 2023
  dplyr::filter(last_yr_st == '2023') 


## Load and process CA lakes

# Load raw CA lake data
CA.Lakes = 
  read_sf(here::here('Data/Geospatial/raw', 
                     'CA_Lakes.shp')) %>%
  # Reproject to NAD_1983_California_Teale_Albers (meters)
  st_transform(., crs = 3310) %>% 
  st_zm(.) %>%                 # remove z dimension
  st_intersection(., CA.bound) # clip lakes to CA boundary
  

# export CA lakes clipped to CA boundary
save(CA.Lakes, 
     file = here::here('Data/Geospatial/', 'CA_Lakes.Rda'))


## Isolate fishing lakes

# find indices of intersecting lakes
Fish.Lakes.idx = 
  st_intersects(Fishing.Locs, CA.Lakes, sparse = TRUE) %>%
  as.numeric()

# Filter to identify fishing lakes

Fish.Lakes = 
  CA.Lakes[Fish.Lakes.idx, ] %>%
  filter(!is.na(LAT_NAD83)) %>%
  dplyr::mutate(
    lake_area_m2 = st_area(.)) %>%  # calculate area of each lake
  dplyr::select(!(REGION:INTPTLON)) # drop UCB cols

# Identify duplicated lakes 

# The CDFW fishing guide provides two point locations for some
# lakes that the CA Lakes shapefile has combined (e.g., Upper
# Echo Lake and Lower Echo lake have separate points in the 
# CDFW fishing guide, but the CA Lakes shapefile combines them 
# into one multiploygon - Upper & Lower Echo Lakes).

Dups.Fish = 
  Fish.Lakes %>%
  st_drop_geometry(.) %>%
  count(DFGWATERID) %>% 
  filter(n > 1) %>%
  pull(DFGWATERID)

# Lines 102-113 is how we determined why lakes were duplicated
# it is only for informational purposes.

#Dup.Lakes = 
#  Fish.Lakes %>% 
#  filter(DFGWATERID %in% Dups.Fish) %>% 
#  mutate(Area = st_area(.))

#check = 
#  Dup.Lakes %>%
#  st_drop_geometry(.) %>%
#  dplyr::select(any_of(c('DFGWATERID', 'NAME', 'LAT_NAD83', 'LON_NAD83'))) %>%
#  unite('coords', 3:4, sep = ', ') 

#write_csv(check, here::here('', 'DupCheck.csv'))

# remove duplicated lakes
Fish.Lakes = 
  distinct(Fish.Lakes, DFGWATERID, .keep_all = TRUE)

# export fishing lakes
save(Fish.Lakes, 
     file = here::here('Data/Geospatial/', 'CA_Fish_Lakes.Rda'))

