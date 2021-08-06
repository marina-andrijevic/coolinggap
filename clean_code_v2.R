
library(ggplot2)
library(tidyverse)
library(broom)
library(countrycode)
library(readstata13)
library(sandwich)
library(lmtest)
library(zoo)
library(plm)
library(broom)
library(readxl)
library(ggthemes)
library(maps)
library(rworldmap)
library(viridis)
library(gridExtra)
library(stargazer)
library(caret)
library(leaps)
library(betareg)
library(mgcv)
library(purrr)


# Modeling steps 


# Decadal Cooling degree days (CDD) data until 2100, for three RCPs (2.6, 4.5 and 6.0) and all SSPs

cdd.18 <- read.csv('data/decadal_cdd/ISO_CDD_GDP18_agg_data.csv') %>% 
  mutate(sptemp = 'cdd18')
cdd.20 <- read.csv('data/decadal_cdd/ISO_CDD_GDP20_agg_data.csv') %>% 
  mutate(sptemp = 'cdd20')
cdd.22 <- read.csv('data/decadal_cdd/ISO_CDD_GDP22_agg_data.csv') %>% 
  mutate(sptemp = 'cdd22')
cdd.24 <- read.csv('data/decadal_cdd/ISO_CDD_GDP24_agg_data.csv') %>% 
  mutate(sptemp = 'cdd24')

cdd <- cdd.18 %>% 
  bind_rows(cdd.20) %>% 
  bind_rows(cdd.22) %>% 
  bind_rows(cdd.24) %>% 
  rename_all(tolower) %>% 
  rename(countrycode = iso) %>%
  filter(stat == 'PopWeightAv' & data == 'CDD') %>% # we will only use the population-weighted CDDs
  rename(cdd = value) %>%
  select(countrycode, sptemp, cdd, rcp, ssp, year) %>% 
  mutate(cdd = coalesce(cdd, 0)) %>% 
  spread(sptemp, cdd) 

cdd.long <- cdd %>% 
  gather(sptemp, cdd, -c(countrycode, rcp, ssp, year)) %>% 
  mutate(sptemp = str_sub(sptemp, 4))


# Take 2020 values to calculate CDD 18 equivalents, to calibrate climate maximum ownership for different set point temperatures
# It would be inaccurate to impose the same climate maximum ownership function on CDDs with set point temperatures above 18oC, 
# i.e. 200 CDDs with the sp temperature of 26oC is not the same as 200 CDDs of 18oC
# In other words, for 200 days of sp 26, you'll need more air conditioning than for 200 days of 18 degrees

cdd.2020 <- cdd %>% 
  filter(year == '2020' & ssp == 'SSP2' & rcp == 'rcp45') %>% 
  select(-c(rcp, ssp, year)) %>% 
  gather(sptemp, cdd, -countrycode) 

cdd.reg <- cdd.2020 %>% spread(sptemp, cdd) 

m1 <- gam(cdd18 ~ s(cdd20), data = cdd.reg) # Fit a spline function to the relationship between CDD 18 and different CDD thresholds
m2 <- gam(cdd18 ~ s(cdd22), data = cdd.reg)
m3 <- gam(cdd18 ~ s(cdd24), data = cdd.reg)

ggplot(cdd.reg, aes(cdd18, cdd24) ) +
  geom_point() +
  stat_smooth(method = gam, formula = y ~ s(x))

# Calculate equivalents for CDDs at other sp temps than 18oC using the functions m1-m3
cdd.eq <- cdd %>% 
  mutate(cddeq20 = predict(m1, newdata = cdd) %>% as.vector(),
         cddeq22 = predict(m2, newdata = cdd) %>% as.vector(),
         cddeq24 = predict(m3, newdata = cdd) %>% as.vector()) %>% 
  mutate(cddeq20 = ifelse(cdd20 == 0, 0, cddeq20),
         cddeq22 = ifelse(cdd22 == 0, 0, cddeq22),
         cddeq24 = ifelse(cdd24 == 0, 0, cddeq24)) %>% 
  as.data.frame() %>% 
  select(countrycode,rcp, ssp, year, cdd18, 9:11) %>% 
  rename(cddeq18 = cdd18) %>% 
  gather(sptemp, cddeq, -countrycode, -rcp, -ssp, -year) %>% 
  mutate(sptemp = str_sub(sptemp, 6, 7)) %>% 
  left_join(cdd.long, by = c('countrycode', 'rcp', 'ssp', 'year', 'sptemp')) %>% 
  mutate(clim.max = 1 - 0.949*exp(-0.00187*cddeq)) 

# Climate maximum ownership as function of CDDs at different set point temperatures

cdd.eq$sptemp <- factor(cdd.eq$sptemp, level  = c('24', '22', '20', '18'))
ggplot(cdd.eq %>% filter(rcp == 'rcp45' & year == 2020 & ssp == 'SSP2')) + 
  geom_line(aes(cdd, clim.max, color = sptemp), size = 1) +
  ylim(0,1) + 
  theme_classic() + 
  scale_color_manual(values = c('18' = '#ffb366', '20'='#ff9933', '22'='#cc6600', '24' = '#804000')) + 
  labs(x = 'Cooling Degree Days', y = 'Climate maximum ownership', color = 'Set-point temperature') +
  theme(legend.text = element_text(size=15),
        legend.position = c(0.8, 0.3),
        legend.background = element_rect(linetype="solid", 
                                         colour ="black"),
        axis.text = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15))
  

# # # # # # # # # # # # # # # # # # # # # # # #
# # # Data on air conditioning coverage
# # # # # # # # # # # # # # # # # # # # # # # #

# Climate maximum ownership based on CDDs of 2010 to be used in regression on observed data

cdd.eq.2010 <- read.csv('Data/cdd_eq_hist.csv') %>% 
  select(-X)

# AC data from different samples; years span 2005-2016

ac.gdl <- read_excel('Data/aircon_gdl.xlsx') %>% 
  mutate(countrycode = countrycode(country, 'country.name', 'iso3c')) %>% 
  group_by(countrycode) %>% 
  filter(year >= 2005 & year <=2016) %>% 
  filter(year == min(year)) %>% # Choose the period and the scenario 
  ungroup() %>% 
  select(countrycode, country, year, airco) 

ac.oecd.epic <- read_excel('Data/ac_oecd_epic.xlsx') # From Randazzo et al. (2020)
  
ac.mastrucci <- read.csv('Data/ac_own_alessio.csv') %>% 
  mutate(country = countrycode(countries, 'iso3c', 'country.name'),
         year = 2010,
         ac_own = ac_own*100) %>%
  rename(countrycode = countries,
         airco = ac_own) %>% 
  filter(!countrycode %in% c('FRA', 'MEX', 'CHN')) # already have the data

ac.iea <- read_excel('Data/Global AC Adoption model_v3-withactuals.xlsx', sheet = 'IEA 2018 data') %>% 
  filter(countrycode %in% c('SAU', 'IDN', 'ZAF')) # Countries we don't have from other sources

ac <- ac.gdl %>% 
  bind_rows(ac.iea) %>% 
  bind_rows(ac.mastrucci) %>%
  bind_rows(ac.oecd.epic) %>% 
  mutate(year = 2010,
         scenario = 'Observed') %>% 
  left_join(cdd.eq.2010, by = 'countrycode') %>% 
  rename(ownership = airco) %>% 
  mutate(avail.ac = ownership/clim.max,
         avail.ac = ifelse(avail.ac > 100, 100, avail.ac)/100)


# AC.availability = AC.ownership/Climate Maximum Saturation

# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Regressions to analyze availability of air conditioning
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Socio-economic covariates 

soc.eco <- read.csv('Data/soc_eco_avg.csv') # Compiled in the ac_data_analysis.R script; variables are average for 2000-2018

master <- ac %>% 
  left_join(soc.eco, by = 'countrycode') %>% 
  select(countrycode, sptemp, clim.max, ownership, avail.ac, avg.gdp, avg.urb, avg.gini, avg.gov) %>% 
  drop_na() %>% 
  mutate(year = 2010) %>% 
  rename(gdppc = avg.gdp,
         gini = avg.gini,
         gov = avg.gov, 
         urbanization = avg.urb)

# Fit a beta regression 
# Group by CDD set point temperatures and fit beta regression models
# We will use different set point temperatures for different countries,
# based on what the lowest residual is in regressions
# Projections will then be done with different set point temperatures too

regions <- read.csv('data/regions.csv')

master.beta <- master  %>% 
  left_join(regions, by = 'countrycode') %>% 
  group_by(sptemp) %>% 
  mutate(avail.ac = (avail.ac*(66)+0.5)/67) %>%  # transform the data to be able to fit a beta regression 
  ungroup() %>% 
  arrange(sptemp)

master.beta.prep <- master.beta %>% 
  filter(!duplicated(countrycode)) 

br1 <- betareg(avail.ac ~ gdppc, data = master.beta.prep, link = 'logit')
br2 <- betareg(avail.ac ~ gdppc + gini, data = master.beta.prep, link = 'logit')
br3 <- betareg(avail.ac ~ gdppc + gini + urbanization, data = master.beta.prep, link = 'logit')

# Tabulated results for Table 1 (SM)
stargazer(br1, br2, br3, 
          title = 'Regression results',
          dep.var.labels = 'Air conditioning availability',
          covariate.labels = c('GDP per capita', 'Inequality', 'Urbanization'),
          omit.stat = "f",
          align = TRUE)


# Fitted values using the regression model br3
countrycode <- master.beta %>% select(countrycode)
fitted.models <- master.beta %>% 
  mutate(id = paste0(countrycode, sptemp)) %>% 
  column_to_rownames(var = 'id') %>% 
  as_tibble() %>% 
  nest(data = -sptemp) %>% # nest the data by its set point temperature and run regressions 
  mutate(fit = purrr::map(data, ~ betareg(avail.ac ~ gdppc + gini + urbanization, data = .x)),
         augmented = purrr::map(fit, augment)) %>% 
  unnest(augmented) %>% 
  bind_cols(countrycode) %>% 
  select(countrycode, sptemp, avail.ac, .fitted)

# Select the country-specific set point temperature for the rest of the analysis based on
# what the smallest residual is 

model.select <- fitted.models %>% 
  group_by(countrycode) %>% 
  mutate(resid = abs(avail.ac - .fitted)) %>% 
  filter(resid == min(resid)) %>% 
  ungroup() %>% 
  arrange(countrycode) %>% 
  mutate(id = paste0(countrycode, sptemp))

# Plot the distribution of sp temperatures that are best fitting 
# Most of them fall into the 24oC threshold, so we will use this estimate for the estiamte of heat stress later on 
ggplot(model.select) + 
  geom_bar(aes(as.factor(sptemp))) +
  labs(x = 'Set point temperature', y = 'Number of countries') +
  theme_classic() + 
  theme(axis.text = element_text(size = 12))

model.id <- model.select$id # These models will be extracted from projections in the end
model.cc <- unique(model.select$countrycode) # Countries from the original sample of ACs


# Validation table for the SM Table 2
validation <- model.select %>% 
  left_join(master %>% select(countrycode, clim.max, sptemp, ownership), by = c('countrycode', 'sptemp')) %>% 
  mutate(ownership.fit = .fitted*clim.max*100, 
         diff = ownership - ownership.fit) %>% 
  mutate(avail.ac = avail.ac*100,
         .fitted = .fitted*100) %>% 
  select(-c(sptemp, id, diff))

#write.csv(validation, 'Manuscript figures/SM_Table2.csv')

# Projections of AC availability based on the coefficient estimates of income and inequality
# and AC ownership based on the formula ac.ownership = climate.max.ownership x ac.availability

cdd.proj <- cdd.eq %>% 
  rename(scenario = ssp)

proj.data <- read.csv('data/soc_eco_projections.csv') %>% 
  right_join(cdd.proj, by = c('countrycode', 'year', 'scenario')) %>% 
  left_join(regions, by = 'countrycode')

yrs <- seq(2010, 2100, 10)

# Projections of ac availability 
proj.fit.prep <- proj.data %>% 
  nest(data = -sptemp) %>% 
  mutate(fit = purrr::map(data, ~betareg(avail.ac ~ gdppc + gini + urbanization, data = master.beta)),
         ac.proj = purrr::map(.x = fit, ~predict(., newdata = proj.data))) %>% 
  select(sptemp, ac.proj) %>% 
  unnest(cols = c(ac.proj)) %>% 
  cbind(proj.data = proj.data) %>% 
  rename(countrycode = proj.data.countrycode, 
         rcp = proj.data.rcp,
         year = proj.data.year, 
         scenario = proj.data.scenario,
         clim.max = proj.data.clim.max) %>% 
  mutate(id = paste0(countrycode, proj.data.sptemp), 
         id2 = paste0(countrycode, sptemp)) 

# Observed sample projections 
proj.fit.1 <- proj.fit.prep %>% 
  filter(id %in% model.id) %>%
  filter(id2 %in% model.id) %>% 
  mutate(ownership.proj = ac.proj*clim.max)

# Projections for all other countries for which we have the socio-economic
# data from the SSPs but not observed AC data

proj.fit.2 <- proj.fit.prep %>% 
  filter(!countrycode %in% model.cc) %>%  # delete countries from the observed sample
  group_by(countrycode, year, rcp, scenario) %>% 
  mutate(ac.proj = mean(ac.proj),
         clim.max = mean(clim.max)) %>% 
  ungroup() %>% 
  mutate(ownership.proj = ac.proj*clim.max) %>% 
  filter(proj.data.sptemp == 24,
         sptemp == 24)

proj.fit <- proj.fit.1 %>% 
  bind_rows(proj.fit.2) %>% 
  select(countrycode, year, rcp, scenario, sptemp, clim.max, ac.proj, ownership.proj)

sample <- sample(proj.fit$countrycode, size = 30)

ggplot() + 
  geom_point(data = ac %>% filter(countrycode %in% sample), aes(year, avail.ac, color = scenario)) + 
  geom_line(data = proj.fit %>% filter(countrycode %in% sample & rcp == 'rcp60'), aes(year, ac.proj*100, color = scenario)) +
  facet_wrap(~countrycode) + 
  labs(x = 'Year', y = 'AC ownership', title = 'RCP6.0', color = 'Scenario') +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12)) +
  scale_color_manual(values = c('SSP1' = '#588c7e', 'SSP2' = '#ffcc5c', 
                               'SSP3' = '#ff6f69', 'SSP4' = '#d96459', 'SSP5' = '#96ceb4', 'Observed' = '#eeac99'))

# Regional aggregation 

ac.regs <- ac %>% 
  left_join(regions, by = 'countrycode') %>% 
  group_by(region) %>% 
  mutate(avg.own = mean(ownership),
         avg.avail = mean(avail.ac, na.rm = T))%>% 
  ungroup() %>% 
  filter(!duplicated(countrycode)) %>% 
  mutate(year = as.integer(year)) %>% 
  drop_na(region) %>% 
  select(countrycode, region, year, scenario, avg.own, avg.avail)
  
proj.regs <- proj.fit %>% 
  left_join(regions, by = 'countrycode') %>% 
  group_by(region, year, scenario) %>% 
  mutate(avg.proj.own = mean(ownership.proj, na.rm = T),
         avg.proj.avail = mean (ac.proj, na.rm = T))%>% 
  ungroup() %>% 
  mutate(year = as.integer(year)) %>% 
  drop_na(ac.proj)


# Regional plot availability & ownership: 

ggplot() + 
  geom_point(data = ac.regs, aes(year, avg.avail*100, color = scenario)) + 
  #geom_point(data = ac.regs, aes(year, avg.own, color = scenario)) +
  geom_line(data = proj.regs %>% filter(rcp == 'rcp45'), aes(year, avg.proj.avail*100, color = scenario)) +
  #geom_line(data = proj.regs %>% filter(rcp == 'rcp45'), aes(year, avg.proj.own*100, color = scenario), linetype = "dashed") +
  facet_wrap(~region, nrow = 2) + 
  labs(x = 'Year', y = 'Share of AC availability in the population', title = 'Regional projections of AC availability (RCP 4.5)', color = 'Scenario') +
  theme(strip.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(fill = "white", colour = "#e5e5e4"),
        legend.key = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = c('SSP1' = '#588c7e', 'SSP2' = '#ffcc5c', 
                                'SSP3' = '#ff6f69', 'SSP4' = '#d96459', 'SSP5' = '#96ceb4', 'Observed' = '#eeac99'))



# Couple with projections of population exposed to heat stress

pop.18 <- read.csv('data/pop_exposure_map/ISO_CDD_pop18_agg_data.csv') %>% 
  mutate(sptemp = 18)
pop.20 <- read.csv('data/pop_exposure_map/ISO_CDD_pop20_agg_data.csv') %>% 
  mutate(sptemp = 20)
pop.22 <- read.csv('data/pop_exposure_map/ISO_CDD_pop22_agg_data.csv') %>% 
  mutate(sptemp = 22)
pop.24 <- read.csv('data/pop_exposure_map/ISO_CDD_pop24_agg_data.csv') %>% 
  mutate(sptemp = 24)

pop.exp <- pop.22 %>% 
  bind_rows(pop.24) %>% 
  bind_rows(pop.20) %>%
  bind_rows(pop.18) %>%
  rename(countrycode = ISO,
         scenario = ssp) %>% 
  mutate(sptemp = as.character(sptemp)) %>% 
  left_join(regions, by = 'countrycode') %>% 
  spread(stat, value) %>% 
  rename(cdd.count = data) %>% 
  group_by(rcp, scenario, year, sptemp, cdd.count) %>% 
  mutate(global.exp = sum(PopSumExposed), 
         global.pop = sum(TotalPop), 
         global.share = global.exp/global.pop) %>% 
  ungroup() %>% 
  filter(!countrycode %in% c('GUF','ESH')) %>% #Don't have the soc-eco data for French Guyana and Western Sahara
  mutate(cdd.count = recode(cdd.count, 'CDD_50' = '>50', 'CDD_100' = '>100', 'CDD_200' = '>200', 'CDD_400' = '>400'))

# Heat map: 
# Heat map with global population (edited in Adobe Illustrator to add representative countries for some of the temperature/CDD count combination)

pop.exp$cdd.count <- factor(pop.exp$cdd.count,levels = c('>50', '>100', '>200', '>400'))
ggplot(pop.exp %>% filter(year == 2050 & scenario == 'SSP2' & rcp == 'rcp45'), aes(y = sptemp, x = cdd.count, fill = global.exp/(10^9))) +
  geom_tile() +
  scale_fill_viridis(option="magma", direction = -1) + 
  labs(x = "Number of CDDs", y = 'Set point temperature', fill = 'Global population\nexposure (billion)\n(2050/SSP2/RCP4.5)') 


# Regional cooling gap based on mean population exposure (from the heat map above)

sptemp.set <- c('24')
cddcount.set <- c('>50', '>100', '>200', '>400')

cool.gap.avg <- pop.exp %>% 
  filter(sptemp %in% sptemp.set & cdd.count %in% cddcount.set) %>%
  group_by(countrycode, rcp, scenario, year) %>% 
  mutate(pop.exp.med = median(PopSumExposed, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(id = paste0(countrycode, year, rcp, scenario)) %>% 
  filter(!duplicated(id)) %>% 
  select(-c(id, sptemp, cdd.count)) %>% 
  right_join(proj.fit %>% select(-sptemp), by = c('countrycode', 'rcp', 'scenario', 'year')) %>% 
  mutate(cooling.gap = pop.exp.med*(1-ac.proj)) %>% 
  group_by(region, rcp, scenario, year) %>% 
  mutate(reg.cool.gap = sum(cooling.gap, na.rm = T),
         reg.pop = sum(TotalPop, na.rm = T),
         reg.share = reg.cool.gap/reg.pop) %>% 
  ungroup() %>% 
  drop_na(region) %>% 
  mutate(rcp = recode(rcp, 'rcp26' = 'RCP2.6', 'rcp45' = 'RCP4.5', 'rcp60' = 'RCP6.0'))

# Bar plots with the average cooling gap
# Only three SSP scenarios but showing the RCP spread
ggplot(cool.gap.avg %>% filter(year == 2100 & scenario %in% c('SSP1', 'SSP2', 'SSP3')), aes(x = region, y = reg.share*100, fill = scenario, color = scenario)) +
  stat_summary(geom = 'bar', fun = mean, position = 'dodge') +
  stat_summary(fun.data = spread.fun, geom = "pointrange", position = position_dodge(width = 0.9), color = '#c9c7c7', size = 0.3) + 
  labs(y = '% of population', x = 'Regions', title = paste0('Share of population affected by cooling gap in 2050'), fill = 'Scenario', color  = 'Scenario') + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 12)) +
  scale_color_manual(values = c('SSP1' = '#efe350ff', 'SSP2' = '#f68f46ff', 
                                'SSP3' = '#a65c85ff')) + 
  scale_fill_manual(values = c('SSP1' = '#efe350ff', 'SSP2' = '#f68f46ff', 
                                'SSP3' = '#a65c85ff'))


# Regional cooling gap with all sptemp/cdd count combinations
cool.gap.all <- pop.exp %>% 
  right_join(proj.fit %>% select(-sptemp), by = c('countrycode', 'rcp', 'scenario', 'year')) %>% 
  mutate(cooling.gap = PopSumExposed*(1-ac.proj)) %>% 
  group_by(region, rcp, scenario, year, sptemp, cdd.count) %>% 
  mutate(reg.cool.gap = sum(cooling.gap, na.rm = T),
         reg.pop = sum(TotalPop, na.rm = T),
         reg.share = reg.cool.gap/reg.pop) %>% 
  ungroup() %>% 
  mutate(rcp = recode(rcp, 'rcp26' = 'RCP2.6', 'rcp45' = 'RCP4.5', 'rcp60' = 'RCP6.0'))


# Global plot & sources of uncertainty

# Central estimate of the cooling gap
central.cool.gap <- cool.gap.avg %>% 
  select(countrycode, year, scenario, rcp, pop.exp.med, ac.proj, cooling.gap, TotalPop) %>% 
  group_by(rcp, scenario, year) %>% 
  mutate(central.cool.gap = sum(cooling.gap, na.rm = T),
         central.pop = sum(TotalPop, na.rm = T),
         central.share = central.cool.gap/central.pop,
         central.pop.exp = sum(pop.exp.med, na.rm = T)) %>% 
  ungroup() %>% 
  select(countrycode, year, scenario, rcp, central.cool.gap, central.pop, central.share, central.pop.exp) %>% 
  filter(rcp == 'RCP4.5' & scenario == 'SSP2' & year %in% c(2050, 2100)) %>% 
  mutate(id = paste0(year, rcp, scenario)) %>% 
  filter(!duplicated(id)) %>%
  mutate(countrycode = 'Global') %>% 
  select(-countrycode, -scenario, -rcp, -id)


global.gap <- cool.gap.all %>% 
  select(countrycode, year, rcp, scenario, cooling.gap, TotalPop, sptemp, cdd.count, PopSumExposed) %>% 
  group_by(rcp, year, scenario, sptemp, cdd.count) %>% 
  mutate(global.cgap = sum(cooling.gap, na.rm = T),
         global.exp = sum(PopSumExposed, na.rm = T),
         global.pop = sum(TotalPop, na.rm = T),
         global.share = global.cgap/global.pop) %>% 
  ungroup() %>% 
  mutate(year = as.integer(year),
         id = paste0(rcp, scenario, year, sptemp, cdd.count)) %>%
  filter(!duplicated(id)) %>% 
  select(-id) %>% 
  mutate(countrycode = 'Global') %>% 
  left_join(central.cool.gap, by = 'year') %>% 
  mutate(diff.gap = global.cgap - central.cool.gap, 
         diff.exp = global.exp - central.pop.exp)

# Uncertainty in cooling gap estimates (deviation from the central estimate)
uncertainty <- global.gap %>% 
  mutate(RCP = ifelse(scenario == 'SSP2' & (sptemp %in% sptemp.set & cdd.count %in% cddcount.set), diff.gap, NA),
         SSP = ifelse(rcp == 'RCP4.5' & (sptemp %in% sptemp.set & cdd.count %in% cddcount.set), diff.gap, NA),
         'Set point temperature' = ifelse(scenario == 'SSP2' & rcp == 'RCP4.5' & (cdd.count %in% cddcount.set), diff.gap, NA),
         'CDD count' = ifelse(scenario == 'SSP2' & rcp == 'RCP4.5' & (sptemp %in% sptemp.set), diff.gap, NA)) %>% 
  select(year, RCP, SSP, 'Set point temperature', 'CDD count') %>% 
  gather(source, estimate, -year) %>% 
  drop_na(estimate) %>% 
  mutate(year = as.factor(year))

uncertainty$source <- factor(uncertainty$source,levels = c('Set point temperature', 'CDD count', 'RCP', 'SSP'))
ggplot(uncertainty %>% filter(year %in% c(2050, 2100))) + 
  geom_boxplot(aes(x = source, y = estimate/10^9, fill = year), color = 'grey50', alpha = 0.8) + 
  labs(x = 'Source of uncertainty', y = 'Population (billion)', fill = 'Year', title = 'Cooling gap uncertainty\n(deviation from the central estimate)') + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=22),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  scale_fill_manual(values = c('2050' = '#077B88', '2100' = '#9DC6CA'))

# Uncertainty in exposed population (deviation from the central estimate)

uncertainty.2 <- global.gap %>% 
  mutate(RCP = ifelse(scenario == 'SSP2' & (sptemp %in% sptemp.set & cdd.count %in% cddcount.set), diff.exp, NA),
         SSP = ifelse(rcp == 'RCP4.5' & (sptemp %in% sptemp.set & cdd.count %in% cddcount.set), diff.exp, NA),
         'Set point temperature' = ifelse(scenario == 'SSP2' & rcp == 'RCP4.5' & (cdd.count %in% cddcount.set), diff.exp, NA),
         'CDD count' = ifelse(scenario == 'SSP2' & rcp == 'RCP4.5' & (sptemp %in% sptemp.set), diff.exp, NA)) %>% 
  select(year, RCP, SSP, 'Set point temperature', 'CDD count') %>% 
  gather(source, estimate, -year) %>% 
  drop_na(estimate) %>% 
  mutate(year = as.factor(year))

uncertainty.2$source <- factor(uncertainty.2$source,levels = c('Set point temperature', 'CDD count', 'RCP', 'SSP'))
ggplot(uncertainty.2 %>% filter(year %in% c(2050, 2100))) + 
  geom_boxplot(aes(x = source, y = estimate/10^9, fill = year), color = 'grey50', alpha = 0.8) +
  labs(x = 'Source of uncertainty', y = 'Population (billion)', fill = 'Year', title = 'Exposed population uncertainty\n(deviation from the central estimate)') + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=22),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  scale_fill_manual(values = c('2050' = '#FB9F8B', '2100' = '#E2CCC6'))



# GDP-CDD-AC plot
scatters <- proj.data %>% 
  select(countrycode,year, scenario, gdppc, rcp, sptemp, cdd) %>% 
  right_join(pop.exp %>% select(countrycode, rcp, scenario, year, cdd.count, sptemp, PopSumExposed, TotalPop), by = c('countrycode', 'year', 'scenario', 'rcp', 'sptemp')) %>% 
  left_join(proj.fit %>% select(countrycode, year, rcp, scenario, ac.proj, ownership.proj), by = c('countrycode', 'year', 'scenario', 'rcp'))

ggplot(scatters %>% filter(year == 2050 & rcp == 'rcp45' & cdd.count == 'CDD_50' & scenario == 'SSP2' & sptemp == '24')) + 
  geom_point(aes(x = log(gdppc), y = ac.proj, size = PopSumExposed/10^9, color = cdd, fill = cdd)) + 
  theme_bw() + 
  scale_fill_viridis(option = 'magma') + 
  scale_color_viridis(option = 'magma') + 
  labs(x = 'GDP per capita', y = "AC availability", size = 'Population exposed\n(billion)', fill = 'CDD', color = 'CDD', 
       title = '2050/SSP2/RCP4.5') + 
  facet_wrap(~sptemp, scales = 'free') 




# Function for box plots 
spread.fun <-function(x){data.frame(ymin=min(x),ymax=max(x),y=mean(x))}

ggplot(cool.gap1 %>% filter(year == 2050 & scenario %in% c('SSP1', 'SSP2', 'SSP3')), aes(x = region, y = reg.share*100, fill = scenario, color = scenario)) +
  stat_summary(geom = 'bar', fun = mean, position = 'dodge') +
  stat_summary(fun.data = spread.fun, geom = "pointrange", position = position_dodge(width = 0.9), color = '#c9c7c7') + 
  labs(y = '% of population', x = 'Regions', title = paste0('2050'), fill = 'Scenario', color  = 'Scenario') + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size=22),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  scale_color_manual(values = c('SSP1' = '#588c7e', 'SSP2' = '#ffcc5c', 
                                'SSP3' = '#ff6f69', 'SSP4' = '#d96459', 'SSP5' = '#96ceb4')) + 
  scale_fill_manual(values = c('SSP1' = '#588c7e', 'SSP2' = '#ffcc5c', 
                               'SSP3' = '#ff6f69', 'SSP4' = '#d96459', 'SSP5' = '#96ceb4'))



 

# As a share of population
plot3 <- function(r, y) {
cool.gap2 %>% 
    filter(rcp == r, year == y) %>% 
    ggplot() +
    geom_boxplot(aes(x = region, y = reg.share, fill = scenario, color = scenario)) +
    labs(y = '% of pop', x = 'Regions', title = paste0('Share of population affected by cooling gap in ', y, ' for ', r)) + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 12)) +
    scale_color_manual(values = c('SSP1' = '#588c7e', 'SSP2' = '#ffcc5c', 
                                'SSP3' = '#ff6f69', 'SSP4' = '#d96459', 'SSP5' = '#96ceb4')) + 
    scale_fill_manual(values = c('SSP1' = '#588c7e', 'SSP2' = '#ffcc5c', 
                               'SSP3' = '#ff6f69', 'SSP4' = '#d96459', 'SSP5' = '#96ceb4'))
}

plot3('RCP4.5', 2100)


# Maps 

# Pop exposed to heat stress in 2020 (median CDD)

years <- c(2050, 2100)
rcps <- c('rcp26', 'rcp60')

for(y in years) {
  for(r in rcps) {
dat <- pop.exp %>% 
  filter(year == y & scenario == 'SSP2' & rcp == r & sptemp == '24' & cdd.count == 'CDD_100') %>% 
  mutate(var = pop.exp.med/TotalPop)

map <- joinCountryData2Map(dat, joinCode = "ISO3", nameJoinColumn = "countrycode")
map_poly <-  fortify(map) %>% 
  merge(map@data, by.x="id", by.y="ADMIN", all.x=T) %>%
  arrange(id, order) %>% 
  mutate(var %>% as.numeric())

temp.map <- ggplot(map_poly, aes( x = long, y = lat, group = group )) +
  coord_map(projection = 'mollweide', xlim = c(-180, 180), ylim = c(-60, 75))  + # Remove antarctica
  geom_polygon(aes(fill = var)) +
  scale_fill_viridis(option = "inferno", direction = -1, na.value="#e1e0e0") +
  labs(fill = 'Value'
       ,title = paste0('Population exposed to heat stress in ', y, ' for ',r)   # Change the title of the map
       ,x = NULL
       ,y = NULL) +
  theme(text = element_text(family = 'Helvetica', color = 'gray40')
        ,plot.title = element_text(size = 18)
        ,axis.ticks = element_blank()
        ,axis.text = element_blank()
        ,axis.line = element_blank()
        ,panel.grid = element_blank()
        ,panel.background = element_rect(fill = 'white')
        ,plot.background = element_rect(fill = 'white')
        ,legend.position = c(.08,.26)
        ,legend.background = element_blank()
        ,legend.key = element_blank()
  ) +
  annotate(geom = 'text'
           ,label = ''
           ,x = 18, y = -55
           ,size = 3
           ,family = 'Helvetica'
           ,color = 'gray50'
           ,hjust = 'left'
  )
ggsave(temp.map, file=paste0("map_", r,"_", y,"_v2.pdf"), width = 16, height = 10, units = "cm")
  }
}


# Figure 1: AC baseline data map

dat <- ac %>% 
  mutate(var = ownership)

map <- joinCountryData2Map(dat, joinCode = "ISO3", nameJoinColumn = "countrycode")
map_poly <-  fortify(map) %>% 
  merge(map@data, by.x="id", by.y="ADMIN", all.x=T) %>%
  arrange(id, order) %>% 
  mutate(var %>% as.numeric())

ggplot(map_poly, aes( x = long, y = lat, group = group )) +
  coord_map(projection = 'mollweide', xlim = c(-180, 180), ylim = c(-60, 75))  + # Remove antarctica
  geom_polygon(aes(fill = var)) +
  scale_fill_viridis(option = "viridis", direction = -1, na.value="#e1e0e0") +
  labs(fill = 'Share of \npopulation'
       ,title = 'AC ownership in the baseline data'   # Change the title of the map
       ,x = NULL
       ,y = NULL) +
  theme(text = element_text(family = 'Helvetica', color = 'gray40')
        ,plot.title = element_text(size = 18)
        ,axis.ticks = element_blank()
        ,axis.text = element_blank()
        ,axis.line = element_blank()
        ,panel.grid = element_blank()
        ,panel.background = element_rect(fill = 'white')
        ,plot.background = element_rect(fill = 'white')
        ,legend.position = c(.08,.26)
        ,legend.background = element_blank()
        ,legend.key = element_blank()
  ) +
  annotate(geom = 'text'
           ,label = ''
           ,x = 18, y = -55
           ,size = 3
           ,family = 'Helvetica'
           ,color = 'gray50'
           ,hjust = 'left'
  )


# Fig 2: Regional projections of AC availability

ggplot() + 
  geom_point(data = ac.regs, aes(year, avg.avail*100, color = scenario)) + 
  geom_line(data = proj.regs %>% filter(rcp == 'rcp45'), aes(year, avg.proj.avail*100, color = scenario)) +
  facet_wrap(~region, nrow = 2) + 
  labs(x = 'Year', y = 'Share of AC availability in the population', title = 'Regional projections of AC availability (RCP 4.5)', color = 'Scenario') +
  theme(strip.background = element_rect(colour = "white", fill = "white"),
        panel.background = element_rect(fill = "white", colour = "#e5e5e4"),
        legend.key = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12), 
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = '#e5e5e4')) +
  scale_color_manual(values = c('SSP1' = '#588c7e', 'SSP2' = '#ffcc5c', 
                                'SSP3' = '#ff6f69', 'SSP4' = '#d96459', 'SSP5' = '#96ceb4', 'Observed' = '#eeac99'))

# Fig 3: Heat stress in 2020, 2050 and 2100 (SSP2/RCP4.5)


# Pop exposed to heat stress in 2020 (median CDD)

years <- c(2020, 2050, 2100)

for(y in years) {
    dat <- cool.gap.avg %>% 
      filter(year == y & scenario == 'SSP2' & rcp == 'RCP4.5') %>% 
      mutate(var = pop.exp.med/TotalPop)
    
    map <- joinCountryData2Map(dat, joinCode = "ISO3", nameJoinColumn = "countrycode")
    map_poly <-  fortify(map) %>% 
      merge(map@data, by.x="id", by.y="ADMIN", all.x=T) %>%
      arrange(id, order) %>% 
      mutate(var %>% as.numeric())
    
    temp.map <- ggplot(map_poly, aes( x = long, y = lat, group = group )) +
      coord_map(projection = 'mollweide', xlim = c(-180, 180), ylim = c(-60, 75))  + # Remove antarctica
      geom_polygon(aes(fill = var)) +
      scale_fill_gradient(low = "yellow", high = "red", na.value = NA) +
      labs(fill = 'Value'
           ,title = paste0('         ', y, '              ')   # Change the title of the map
           ,x = NULL
           ,y = NULL) +
      theme(text = element_text(family = 'Helvetica', color = 'gray40')
            ,plot.title = element_text(size = 18)
            ,axis.ticks = element_blank()
            ,axis.text = element_blank()
            ,axis.line = element_blank()
            ,panel.grid = element_blank()
            ,panel.background = element_rect(fill = 'white')
            ,plot.background = element_rect(fill = 'white')
            ,legend.position = c(.08,.26)
            ,legend.background = element_blank()
            ,legend.key = element_blank()
      ) +
      annotate(geom = 'text'
               ,label = ''
               ,x = 18, y = -55
               ,size = 3
               ,family = 'Helvetica'
               ,color = 'gray50'
               ,hjust = 'left'
      )
    ggsave(temp.map, file=paste0("fig2_rcp45","_", y,"heat.pdf"), width = 16, height = 10, units = "cm")
  }


# Fig X Regional cooling gap 

# Bar plots with the average cooling gap
# Only three SSP scenarios but showing the RCP spread
ggplot(cool.gap.avg %>% filter(year == 2100 & scenario %in% c('SSP1', 'SSP2', 'SSP3')), aes(x = region, y = reg.share*100, fill = scenario, color = scenario)) +
  stat_summary(geom = 'bar', fun = mean, position = 'dodge') +
  stat_summary(fun.data = spread.fun, geom = "pointrange", position = position_dodge(width = 0.9), color = '#c9c7c7', size = 0.3) + 
  labs(y = 'Population share (%)', x = 'Regions', title = paste0('2100'), fill = 'Scenario', color  = 'Scenario') + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(color = '#e5e5e4'),
        panel.grid.major.y = element_line(color = '#e5e5e4')) +
  scale_color_manual(values = c('SSP1' = '#F6A97A', 'SSP2' = '#D44292', 
                                'SSP3' = '#4B2991')) + 
  scale_fill_manual(values = c('SSP1' = '#F6A97A', 'SSP2' = '#D44292', 
                               'SSP3' = '#4B2991'))



# Fig X Global figure 

# Global plot
global.plot <- cool.gap.avg %>% 
  select(countrycode, year, scenario, rcp, pop.exp.med, ac.proj, cooling.gap, TotalPop) %>% 
  group_by(rcp, scenario, year) %>% 
  mutate(global.cool.gap = sum(cooling.gap, na.rm = T),
         global.pop = sum(TotalPop, na.rm = T),
         global.share = global.cool.gap/global.pop,
         global.pop.exp = sum(pop.exp.med, na.rm = T)) %>% 
  ungroup() %>% 
  select(countrycode, year, scenario, rcp, global.cool.gap, global.pop, global.share, global.pop.exp) %>% 
  mutate(id = paste0(year, rcp, scenario)) %>% 
  filter(!duplicated(id)) %>%
  mutate(countrycode = 'Global')


ggplot(global.plot %>% filter(rcp == 'RCP4.5')) +
  geom_boxplot(aes(x = as.factor(year), y = global.cool.gap/10^9), fill = '#7c3494', color = 'grey50', alpha = 0.8) +
  labs(y = 'Population (billion)', x = 'Year', title = 'Global population affected by cooling gap in all SSPs (RCP 4.5)') + 
  theme_bw() + 
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

ggplot(global.plot %>% filter(rcp == 'RCP4.5')) +
  geom_boxplot(aes(x = as.factor(year), y = global.pop/10^9))


sm.table1 <- regions %>% 
  right_join(ac %>% select(countrycode), by = 'countrycode') %>% 
  filter(!duplicated(countrycode)) %>% 
  mutate(country = countrycode(countrycode, 'iso3c', 'country.name')) %>% 
  arrange(region) %>% 
  select(country, region) %>% 
  drop_na(country)

write.csv(sm.table1, 'Manuscript figures/smtable.csv')


# Temp check 

tw.days <- read.delim('data/gtc_24.6_days_EWEMBI_2000_2013.txt') %>% 
  mutate(days = parse_number(country)) %>% 
  mutate(countrycode = countrycode(country, 'country.name', 'iso3c')) %>% 
  select(countrycode, days) %>% 
  filter(!days >= 10^20)


tw.merge <- read.csv('Data/tw_mean_allcntry.csv') %>% 
  mutate(countrycode = countrycode(ISO, 'iso2c', 'iso3c')) %>% 
  right_join(cdd.eq.2010, by = 'countrycode') %>% 
  right_join(tw.days, by = 'countrycode') %>% 
  drop_na(sptemp)

ggplot(tw.merge %>% filter(sptemp == 24 & cdd > 50)) +
  geom_point(aes(x = days, y = cdd)) + 
  facet_wrap(~sptemp) + 
  labs(x = 'No of days above 24.6', y = 'CDD')


# Pop size covered by the sample

pop.size <- read.csv('/Users/marinaandrijevic/Documents/SSP data/pop_size_obs_proj.csv', skip = 8) %>% 
  filter(Year == 2010 & Scenario == 'SSP2') %>%
  mutate(countrycode = countrycode(Area, 'country.name', 'iso3c')) %>% 
  filter(!Area == 'World') %>% 
  mutate(pop.total = sum(Population)) %>% 
  right_join(master.beta.prep %>% filter(sptemp == 18), by = 'countrycode') %>% 
  drop_na(countrycode) %>% 
  mutate(sample.pop = sum(Population),
         share = sample.pop/pop.total) 

head(pop.size %>% filter(Area == 'World'))

6972709