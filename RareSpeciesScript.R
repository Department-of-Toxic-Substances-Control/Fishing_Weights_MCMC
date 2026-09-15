#load packages and set working drive
library(sf)
library(tidyverse)
library(cowplot)
library(ggspatial)
library(ggpubr)

here::i_am('Fishing_Weights_MCMC.Rproj')

##############################################################
# This script is used to examine the co-location of popular  #
# fishing lakes with areas of relatively high rare species   #
# richness. It is also used to produce Fig 1 in the main     #
# text.                                                      #
##############################################################

## Load data

# import CDFW Terrestrial species of special concern
rarerich = 
  read_sf(here::here('Data/Geospatial/raw/',
                     'Terrestrial_Biodiversity_Summary_-_ACE_[ds2739].shp')) %>%
  # Reproject to NAD_1983_California_Teale_Albers (meters)
  st_transform(., crs = 3310)

# create subset containing only hexagons with Rare Species ecoregion rank of 4 or higher
RarRankEco45 = subset(rarerich, RarRankEco >= 4)

# load fishing lakes
load(here::here('Data/Geospatial/', 'CA_Fish_Lakes.Rda'))

# load UCB CA state boundary
load(here::here('Data/Geospatial/', 'tl_2023_California.Rda'))

# intersect RarRankEco45 with Fishing lakes file
intersects = st_intersection(RarRankEco45, Fish.Lakes)

## Generate counts

# Calculate percentage of RarrankEco >= 4  that contain or touch a fishing lake. 
n_distinct(intersects$Hex_ID) # answer: 618
n_distinct(RarRankEco45)      # answer 21572
618/21572                     # 0.0286 or 3%

# calculate number of fishing lakes
n_distinct(Fish.Lakes$DFGWATERID) # answer 278

# calculate number of fishing lakes that intersect with RarRankEco >= 4
n_distinct(intersects$DFGWATERID) # answer 200 


## Create rare species map
rare.species.map = 
  ggplot() +
  geom_sf(data = CA.bound, fill = '#ffffff') +
  geom_sf(data = RarRankEco45, aes(fill = as.factor(RarRankEco)), 
          colour = NA) +
  geom_sf(data = st_centroid(Fish.Lakes), aes(colour = 'Fishing Lake'),  
          size = 0.35, show.legend = 'point') +
  scale_fill_manual(
    name = 'Normalized Rare\nSpecies Richness Rank',
    values = c('4' = '#a1d99b', '5' = '#31a354'),
    guide = 
      guide_legend(
        override.aes = list(linetype = "blank", shape = NA))) +
  scale_colour_manual(
    values = c('Fishing Lake' = '#252525'), name = NULL,
    guide = 
      guide_legend(
        override.aes = list(linetype = "blank", shape = 16, size = 3))) +
  annotation_scale(location = 'bl') +
  theme_classic2()

rare.species.map

## Make inset map
NA.countries = 
  read_sf(here::here('Data/Geospatial/raw/', 'boundaries_p_2021_v3.shp'))
    
states = c('US-AZ', 'US-CA', 'US-ID', 'US-OR', 'US-NV', 'US-UT')  
  
display_window = 
  st_sfc(st_point(c(-125, 31)), st_point(c(-113, 43)),
         crs = 4326) %>%
  st_coordinates(.)

MX = 
  NA.countries %>% 
  filter(COUNTRY == 'MEX') %>%
  st_union(.) %>%
  st_transform(., crs = 4326)

W.USA = 
  NA.countries %>%
  filter(STATEABB %in% states) %>%
  mutate(col = if_else(STATEABB == 'US-CA', '1', '0')) %>%
  st_transform(., crs = 4326)
  

# make tibble of major cities
cities = 
  tibble('name' = c('San Francisco', 'Los Angeles', 'San Diego'),
         'x' = c(-122.416389,-118.25, -117.1625), 
         'y' = c(37.7775, 34.05, 32.715)) %>% 
  st_as_sf(., coords = c('x', 'y'), crs = 4326)

NA.inset.map = 
  ggplot() +
  geom_sf(data = W.USA, aes(fill = col), colour = '#4d4d4d') +
  geom_sf(data = MX, fill = '#dfc27d', colour = '#4d4d4d') +
  geom_sf(data = cities, colour = '#b2182b', size = 1.5) +
  geom_sf_text(data = cities, aes(label = name), 
               size = 3.5, colour = 'black', angle = 3,
               nudge_x = c(2, 1.5, 1.5), 
               nudge_y = c(0.5, 0.5, 0.5)) +
  annotate(geom = 'text', x = -123.5, y = 35.5, label = 'Pacific\nOcean') +
  annotate(geom = 'text', x = -115.5, y = 32.25, label = 'Mexico', angle = 3) +
  scale_fill_manual(values = c('#e0e0e0', '#ffffff'), guide = 'none') +
  coord_sf(xlim = display_window[,'X'], 
          ylim = display_window[,'Y'],
         datum = 4326, expand = FALSE) +
  theme_classic2() +
  theme(
    panel.background = element_rect(fill = '#4393c3'),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 1))

species.legend = get_legend(rare.species.map)

## Combine all plots
Fig1 = 
  ggdraw(rare.species.map + theme(legend.position = 'none')) + 
  draw_plot(NA.inset.map, x = 0.30, y = 0.22, scale = 0.35) + 
  draw_plot(species.legend, x = -0.3, y = -0.3, scale = 1.1)

ggsave(here::here('Figures/', 'Fig1.png'), plot = Fig1,
       width = 8.5, height = 11, units = 'in', dpi =300)
