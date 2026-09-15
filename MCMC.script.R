library(readxl)
library(fitdistrplus)
library(ggpubr)
library(bestNormalize)
library(sf)
library(wiqid)
library(magrittr)
library(here)  
library(reshape2)
library(tidyverse)

here::i_am('Fishing_Weights_MCMC.Rproj')


##############################################################
# This script builds and runs the Markov Chain Monte Carlo   #
# simulation.gut check on the MCMC estimates. #
# The exported table is Table S3 in the ESM.                 #
##############################################################

## Statistics

### Import Sinker Loss Data

# read in weight loss data
PbWeightLoss.dat.raw = 
  read_excel(here::here('Data/', 'LFW_loss_dat.xlsx'),
             sheet = 'AllData',
             col_types = c('text', 'text', 'text',
                           'text', 'text', 'numeric',
                           'numeric', 'numeric', 'text'))

### Determine Data Distribution

# visualize distribution of the data
plotdist(PbWeightLoss.dat.raw$Sinkers_nm2, histo = TRUE, demp = TRUE)

#### Transform Sinker Data

# create input vector of sinker loss
sinker.dat = 
  PbWeightLoss.dat.raw %>%
  dplyr::mutate(
    across(Sinkers_nm2, ~if_else(.x == 0, 0.000001, .x))) %>%
  dplyr::pull(Sinkers_nm2)

# run bestNormalize to transform data to norm dist
BN.sinker = bestNormalize(sinker.dat)

# verify transformation worked

BNx = predict(BN.sinker) # transform data

BNx2 = predict(BN.sinker, newdata = BNx, 
               inverse = TRUE) # reverse transformation

all.equal(BNx2, sinker.dat) # verify transformation 


### Markov Chain Monte Carlo Simulation

# run MCMC on transformed data
MCMC.sinker = 
  Bnormal(BN.sinker$x.t,
          #priors=NULL, 
          chains = 10,
          draws = 1e6,           # 1,000,000 draws
          burnin = 5000,         # taken from https://doi.org/10.1021/es304331m
          thin = 15,             # taken from https://doi.org/10.1021/es304331m
          adapt = 1000, 
          doPriorsOnly = FALSE,  # return posterior distributions
          parallel = TRUE,
          seed = 2013)           # set a seed for reproducibility


#### MCMC diagnostics

mcmcOutput::diagPlot(MCMC.sinker, params = 'mu')

# The code below replicates the output of diagPlot above in ggplot
# to create print quality figures. All code was taken from the 
# helper functions in the wiqid package.

nChains = attr(MCMC.sinker, 'nChains')

mat = matrix(MCMC.sinker[,1], ncol = nChains)
name = colnames(MCMC.sinker)[1]

##### Create print quality traceplot
dat = as.data.frame(mat) 
dat$id = 1:nrow(dat)

dat.plt = melt(dat, id.var = 'id')

trace.plt = 
  ggplot(dat.plt, aes(x = id, y = value, 
                    group = variable, colour = variable)) +
  geom_line(alpha = 0.7) +
  scale_colour_brewer(palette = 'Spectral', guide = 'none') +
  geom_hline(yintercept = mean(mat), linetype = 2) + 
  labs(title = 'mu: Rhat = 1', # as calculated by diagPlot()
       x = 'Iterations', 
       y = 'Transformed mu') + 
  theme_classic2() +
  theme(axis.title = element_text(size = 14),
        plot.title = element_text(face = 'bold')) 

  
##### Create print quality density plot

bw = bw.nrd0(mat)

# helper fnc to fit folded density from mcmcOutput
densFold0 = function(mat, bw, from=NA, to=NA,...)  {
  
  # deal with folding for probability and non-negative values
  # use these values if folding is not needed:
  if(is.na(from))
    from = min(mat) - 3*bw
  if(is.na(to))
    to = max(mat) + 3*bw
  mult = 1
  xx = mat
  
  if (min(mat) >= 0 && min(mat) < 2 * bw) {  # it's non-negative
    from = 0
    xx = rbind(mat, -mat)
    mult = 2
  }
  if (min(mat) >= 0 && max(mat) <= 1 &&
      (min(mat) < 2 * bw || 1 - max(mat) < 2 * bw)) { # it's a probability
    to = min(to, 1)
    xx = rbind(mat, -mat, 2-mat)
    mult = 3
  }
  
  # fit density to each column
  n = 512
  dens = apply(xx, 2, function(x) density(x, bw=bw, from=from, to=to, n=n)$y)
  
  return(list(
    x = seq(from, to, length.out=n),
    y = dens * mult))
}

dens = densFold0(mat = mat, bw = bw, from = NA, to = NA)
 
dens.dat = as.data.frame(dens)
dens.dat$id = 1:nrow(dens.dat)


dens.dat.plt = 
  dens.dat %>%
  pivot_longer(
    cols = starts_with('y.'),
    names_to = 'variable',
    values_to = 'value')

dens.plt = 
  dens.dat.plt %>%
  ggplot() +
  geom_line(aes(x = x, y = value, 
                group = variable, colour = variable), alpha = 0.7) +
  scale_colour_brewer(palette = 'Spectral', guide = 'none') +
  geom_hline(yintercept = 0, colour = '#cccccc') + 
  geom_vline(xintercept = mean(dens.dat.plt$x), linetype = 2, colour = '#de2d26') + 
  labs(title = 'MCEpc = 0.03', # as calculated by diagPlot()
         x = 'Transformed mu', 
         y = 'Density') + 
  theme_classic2() +
  theme(axis.title = element_text(size = 14),
        plot.title = element_text(face = 'bold')) + 
  scale_x_continuous(breaks = c(-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6))

##### create combined figure for ESM

diag.plot = 
  ggarrange(trace.plt, dens.plt, ncol = 2, widths = c(2, 1),
            labels = c('a.', 'b.'), label.x = c(0.15, 0.22), label.y = 0.9)

ggsave(here::here('Figures', 'FigS1.png'), plot = diag.plot,
       width = 6.5, height = 3, units = 'in', dpi =300)

#### MCMC Estimates

# extract high, median, and low scenarios
MCMC.out = 
  tibble(mu_trans = MCMC.sinker$mu, 
         s_trans = MCMC.sinker$sigma) %>% 
  dplyr::mutate(
    mu = predict(BN.sinker, newdata = mu_trans, inverse = TRUE),
    s = predict(BN.sinker, newdata = s_trans, inverse = TRUE))

loss.rates = 
  quantile(MCMC.out$mu, c(0.1, 0.5, 0.9)) %>% 
  tibble::enframe(name = 'scenario', value = 'sinker_dens_u') %>% 
  left_join( 
    quantile(MCMC.out$s, c(0.1, 0.5, 0.9)) %>% 
      tibble::enframe(name = 'scenario', value = 'sinker_dens_sd'),
    by = c('scenario' = 'scenario'))


## Calculate Statewide Losses

### Load Geospatial Data

# load Californa fishing lakes
load(here::here('Data/Geospatial/', 'CA_Fish_Lakes.Rda'))

### Load Retailer Data
retailer.dat = 
  read_csv(here::here('Data/', 'RetailerSinkers.csv'))

# generate summary stats for retailer data
retail.sd = signif(sd(retailer.dat$g), digits = 3)
retail.med = signif(median(retailer.dat$g), digits = 3)

# generate table S2
Retail.Summ =
  retailer.dat %>%
  group_by(Size) %>% 
  summarise(mass = signif(min(g), digits = 3),
            ct = n()) %>%
  bind_rows(
    tibble(Size = 'Median', mass = retail.med),
    tibble(Size = 'Std. Dev.', mass = retail.sd))

# export table
write_csv(Retail.Summ, here::here('Figures/', 'TableS2.csv'))

### Calc losses
Sinker.Lakes = 
  Fish.Lakes %>%
  dplyr::mutate(
    low_scenario = 
      signif(lake_area_m2 * as.numeric(loss.rates[1,2]), 
             digits = 3),
    med_scenario = 
      signif(lake_area_m2 * as.numeric(loss.rates[2,2]),
             digits = 3),
    high_scenario = 
      signif(lake_area_m2 * as.numeric(loss.rates[3,2]),
             digits = 3)) %>%
  relocate(COUNTY, .after = NAME) %>% # move county to make cleaning easier
  dplyr::select(!(5:15))

County.Summ = 
  Sinker.Lakes %>% 
  st_drop_geometry() %>%
  dplyr::mutate(
    across(ends_with('_scenario'), ~units::drop_units(.x))) %>%
  group_by(COUNTY) %>%
  dplyr::summarise(
    low_scenario_total = 
      signif(sum(low_scenario, na.rm = TRUE), digits = 3),
    med_scenario_total = 
      signif(sum(med_scenario, na.rm = TRUE), digits = 3),
    high_scenario_total = 
      signif(sum(high_scenario, na.rm = TRUE), digits = 3)) %>%
  bind_rows(summarise(., 
                      across(where(is.numeric), sum),
                      across(where(is.character), ~'Total')))

# generate data for table 1
Loss.Summ = 
  County.Summ %>%
  slice_tail(n = 1) %>% 
  pivot_longer(cols = ends_with('total'),
               names_to = 'Scenario',
               values_to = 'Number of Weights') %>% 
  dplyr::select(-COUNTY) %>%
  mutate(
    across(1, ~str_remove_all(.x, '(?<=_).*')),
    across(1, ~str_remove_all(.x, '_')),
    `Mass of Pb` = 
      signif(`Number of Weights` * 7.09e-6, digits = 3))

# export table
write_csv(Loss.Summ, here::here('Figures/', 'Table1.csv'))

### Sensitivity Analysis
Sens.Anal = 
  Loss.Summ %>%
  mutate(
    `Mass of Pb low` = 
      signif(`Number of Weights` * ((retail.med - retail.sd)*1e-6),
             digits = 3),
    `Mass of Pb hgh` = 
      signif(`Number of Weights` * ((retail.med + retail.sd)*1e-6), 
             digits = 3)) %>%
  dplyr::select(-2) %>%
  pivot_longer(
    cols = starts_with('Mass'),
    names_to = 'Type',
    values_to = 'Mass') %>%
  mutate(across(Type, 
                ~case_when(
                  str_detect(.x, 'low') ~ 'median - 1 SD',
                  str_detect(.x, 'hgh') ~ 'median + 1 SD',
                  .default = 'median')),
         across(Scenario, 
                ~factor(.x, levels = c('low', 'med', 'high'))),
         across(Type, 
                ~factor(.x, 
                        levels = c('median - 1 SD', 
                                   'median', 
                                   'median + 1 SD')))) %>%
  arrange(Type, Scenario) %>%
  mutate(Year_Loss = 
           signif(Mass/73, digits = 3),
         Angler_Loss = 
           signif((Year_Loss*1000)/1440000, digits = 3)) #1,440,000 is the number of freshwater anglers

Sens.Anal.plt.dat = 
  Sens.Anal %>% 
  dplyr::select(c(1:3)) %>%
  pivot_wider(names_from = Type,
              values_from = Mass) %>%
  setNames(c('Scenario', 'lower', 'Mass', 'upper'))
  

# visualize sensitivity analysis
Sens.plot.a = 
  ggplot() +
  geom_errorbar(data = Sens.Anal.plt.dat,
                mapping = aes(xmin = lower, xmax = upper, y = Scenario), 
                width = 0.15, linewidth = 0.5, linetype = 5, 
                colour = '#737373') +
  geom_point(data = Sens.Anal.plt.dat, 
             aes(x = Mass, y = Scenario, fill = Scenario), shape = 21, 
             size = 2.5, stroke = 0.75) +
  scale_fill_brewer(palette = 'RdBu', direction = -1, 
                    guide = 'none') +
  labs(title = '', 
       x = 'Mass (tonnes Pb)', 
       y = 'MCMC Scenario') + 
  theme_classic2() +
  theme(axis.title = element_text(size = 12),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        plot.margin = margin(5.5, 2.5, 5.5, 1.5))

# make a df of per angler losses for comparison to other studies
Compares = 
  #tibble(Estimate 
  #       = c('U.S. annual losses', 'Market sales', 
  #           'EU25 – low', 'EU25 – median', 'EU25 – high'),
  #       Angler_Loss = c(0.114, 0.142, 0.1, 0.2, 0.301)) %>%
  tibble(Estimate 
         = c('U.S. annual losses', 'Market sales', 'EU25'),
         Angler_Loss = c(0.114, 0.142, 0.2),
         lower = c(NA, NA, 0.1),
         upper = c(NA, NA, 0.301)) %>%
  bind_rows(Sens.Anal %>% 
              dplyr::select(-c(3:4)) %>%
              pivot_wider(names_from = Type, 
                          values_from = Angler_Loss) %>%
              #unite('Estimate', Scenario:Type, remove = TRUE) %>%
              setNames(c('Estimate', 'lower', 'Angler_Loss', 'upper')) %>%
              mutate(across(Estimate, 
                            ~str_c('MCMC – ', .x)))) %>%
  mutate(
    #across(Estimate, 
    #       ~factor(.x, 
    #               levels = c('MCMC – low', 'MCMC – med', 'MCMC – high',
    #                          'Market sales', 'U.S. annual losses',
    #                          'EU25 – low', 'EU25 – median', 'EU25 – high'))),
    across(Estimate, 
           ~factor(.x, 
                   levels = c('MCMC – low', 'MCMC – med', 'MCMC – high',
                              'Market sales', 'U.S. annual losses',
                              'EU25'))),
    lab = case_when(
      str_detect(Estimate, 'MCMC') ~ 'MCMC',
      str_detect(Estimate, 'EU')   ~ 'EU',
      .default = as.character(Estimate)))

Sens.plot.b = 
  Compares %>%
  ggplot(aes(x = Angler_Loss, y = Estimate)) +
  geom_errorbar(data = Compares,
                mapping = aes(xmin = lower, xmax = upper, y = Estimate), 
                width = 0.15, linewidth = 0.5, colour = '#737373', 
                linetype = 2) +
  geom_point(aes(fill = lab), shape = 21, size = 2.5, stroke = 0.75) +
  scale_fill_brewer(palette = 'Pastel1', guide = 'none') +
  labs(title = '', 
       x = 'Per Angler Loss (kg Pb/yr)', 
       y = '') + 
  theme_classic2() +
  theme(axis.title = element_text(size = 12),
        plot.margin = margin(5.5, 5.5, 5.5, .25))

# Create combined plot for SI
Sens.plot = 
  ggarrange(Sens.plot.a, Sens.plot.b, ncol = 2, widths = c(1.25, 1.75),
            labels = c('a.', 'b.'), 
            label.x = c(0.21, 0.37), label.y = 0.9)

  
ggsave(here::here('Figures', 'FigS2.png'), plot = Sens.plot,
       width = 6.5, height = 3, units = 'in', dpi =300)


# export table
write_csv(Sens.Anal %>%
            setNames(c('Scenario', 'Type', 'Total Mass (tonnes)', 
                       'Total Yearly Loss (tonnes Pb/year)', 
                       'Per Angler Loss (kg Pb/year)')), 
          here::here('Figures/', 'TableS3.csv'))


