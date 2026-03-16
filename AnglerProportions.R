library(tidyverse)

here::i_am('Fishing_Weights_MCMC.Rproj')

##############################################################
# This script is to run the gut check on the MCMC estimates. #
# The exported table is Table S3 in the ESM.                 #
##############################################################

## Load data

# CA Fishing License sales from 2010 - 2019 (CDFW 2021)
# CA Fishing License sales from 2020 - 2022 (CDFW 2023b)
CDFW.dat = 
  tibble(
    'Licenses' = c(1811311, 1742782, 1803632, 1811547, 1782060, 
                   1766030, 1790870, 1788354, 1776844, 1643348, 
                   1959119, 1840018, 1642994),
    'Lifetimes' = c(942, 987, 1078, 1195, 1195, 1276, 1341, 
                    1500, 1753, 1772, 2302, 2443, 2374))

# U.S. Fishing Licenses (millions) from 2010 - 2021 (Statista 2022)
# U.S. anglers (millions) 2022 (U.S. FWS 2022)
US.dat = 
  tibble(
    'US_anglers' = c(36.3, 36.3, 37.6, 36.0, 37.1, 38.1, 38.6, 
                     38.4, 39.2, 41.4, 38.5, 39.3, 40.0))

## Calculate Proportion of CA:US angles 
Prop.calcs = 
  tibble('Year' = seq.int(2010, 2022, by =1)) %>%
  bind_cols(CDFW.dat, US.dat) %>%
  # Calc total CA anglers (millions)
  mutate(    
    CA_anglers = signif((Licenses + Lifetimes)/10^6, digits = 3), 
    .after = Year) %>%
  # Calc proportion of CA anglers to U.S. Anglers
  mutate(    
    Prop_CA = signif(CA_anglers/US_anglers, digits = 3)) %>%
  dplyr::select(-c(3:4))

write_csv(Prop.calcs, here::here('Figures/', 'TableS3.csv'))
