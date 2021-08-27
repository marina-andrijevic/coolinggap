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


# Observed data for estimating the regression model
# Averaged for 2005-2016 to smooth possible yearly volatilities when matching the data with the AC data

gdp <- read.csv("Data/gdppc_wb.csv", skip = 3) %>%
  rename_all(tolower) %>%
  pivot_longer(x1960:x2019, names_to = "year", values_to = "gdppc") %>%
  mutate(year = year %>% str_replace("x", "") %>% as.integer()) %>%
  select(country.code, year, gdppc) %>%
  filter(year >= 2005 & year <= 2016) %>%
  rename(countrycode = country.code) %>%
  group_by(countrycode) %>%
  mutate(avg.gdp = mean(gdppc, na.rm = T)) %>%
  ungroup() %>%
  filter(!duplicated(countrycode)) %>%
  select(-year)

governance <- read.csv('Data/governance_readiness_hist_proj.csv') %>% 
  filter(year >= 2005 & year <= 2016) %>%
  select(-X, -scenario) %>% 
  group_by(countrycode) %>% 
  mutate(avg.gov = mean(governance, na.rm = T)) %>% 
  ungroup() %>% 
  filter(!duplicated(countrycode)) %>% 
  select(-year)

urb <- read_excel('Data/WUP2018-F02-Proportion_Urban.xls', skip = 16) %>% 
  rename_all(tolower) %>% 
  pivot_longer('2005':'2015', names_to = "year", values_to = "urban") %>% 
  select(countrycode, year, urban) %>% 
  filter(countrycode < 900) %>% 
  mutate(countrycode = countrycode(countrycode, 'iso3n', 'iso3c')) %>% 
  group_by(countrycode) %>% 
  mutate(avg.urb = mean(urban, na.rm = T)) %>% 
  ungroup() %>%
  filter(!duplicated(countrycode)) %>% 
  select(-year)

ineq <- read_excel('Data/WIID_inequality.xlsx') %>% 
  select(c3, year, gini_reported) %>% 
  filter(year >= 2005 & year <= 2016) %>% 
  rename(countrycode = c3) %>% 
  group_by(countrycode) %>% 
  mutate(avg.gini = mean(gini_reported, na.rm = T)) %>% 
  ungroup() %>% 
  filter(!duplicated(countrycode)) %>% 
  select(-year)

iam.regs <- read_excel('Data/iam_regions.xlsx') %>% 
  mutate(countrycode = countrycode(country, 'country.name', 'iso3c'))
                   
soc.eco <- gdp %>%  
  left_join(governance, by = 'countrycode') %>% 
  left_join(urb, by = 'countrycode') %>% 
  left_join(ineq, by = 'countrycode')

write.csv(soc.eco, 'Data/soc_eco_avg.csv')


# Projections data 

scen <- data.frame(scenario = c('SSP1', 'SSP2', 'SSP3', 'SSP4', 'SSP5'))

gdp <- read.csv('Data/gdp_yearly.csv')

gdp1020 <- gdp %>% 
  filter(year == 2020) %>% 
  select(countrycode, year, gdppc) %>% 
  merge(scen)

gini1020 <- read_excel('data/Gini_projections_SSPs.xlsx', sheet = 'projected_ginis_full-set') %>% 
  gather(countrycode, gini, -scenario, -year) %>% 
  filter(year == 2020 & scenario == 'SSP2') %>% 
  select(countrycode, year, gini) %>% 
  merge(scen)

gini.proj <- read_excel('data/Gini_projections_SSPs.xlsx', sheet = 'projected_ginis_full-set') %>% 
  gather(countrycode, gini, -scenario, -year) %>% 
  filter(year > 2020) %>% 
  bind_rows(gini1020)

urb.1020 <- read.csv('data/urb_projections.csv', sep = ';') %>% 
  rename(scenario = scenaio) %>% 
  gather(year, urbanization, -scenario, -countrycode) %>% 
  mutate(year = year %>% str_replace("X", "") %>% as.integer()) %>%
  filter(year == 2020 & scenario == 'SSP2') %>% 
  select(countrycode, year, urbanization) %>% 
  merge(scen)

urb.proj <- read.csv('data/urb_projections.csv', sep = ';') %>% 
  rename(scenario = scenaio) %>% 
  gather(year, urbanization, -scenario, -countrycode) %>% 
  mutate(year = year %>% str_replace("X", "") %>% as.integer()) %>% 
  filter(year > 2020) %>% 
  bind_rows(urb.1020)

soc.eco.proj <- read.csv('data/gdp_edu_projections.csv') %>% 
  select(countrycode, year, scenario, gdppc) %>% 
  filter(year > 2020) %>% 
  bind_rows(gdp1020) %>% 
  left_join(gini.proj, by = c('countrycode', 'year', 'scenario')) %>%
  left_join(urb.proj, by = c('countrycode', 'year', 'scenario'))

write.csv(soc.eco.proj, 'Data/soc_eco_projections.csv')

