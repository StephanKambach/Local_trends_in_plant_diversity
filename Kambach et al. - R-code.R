########################################
# Analyse local changes in plant diversity from ReSurveyEurope Data
# Script by: Stephan Kambach using code from Helge Bruelheide
# stephan.kambach@gmail.com, helge.bruelheide@botanik.uni-halle.de
# initial date: 28.06.2025
#######################################

# clear workspace
rm(list=ls())
gc()

# load libraries
library(broom)
library(broom.mixed)
library(codyn)
library(data.table)
library(DescTools)
library(FD)
library(lme4)
library(lubridate)
library(modEvA)
library(readxl)
library(tidymodels)
library(tidyverse)
library(usedist)
library(vegan)
library(weights)
library(picante)

# set working directory
setwd("")

# load additional functions
source("Kambach et al. - additional functions.R")

# 1. load data ------------------------------------------------------------

# load ReSurveyEurope data
dat.resurvey = readRDS("dat_resurvey.RDS")

# load hand-assigned habitat-change trajectories 
# (stable, succession, disturbance, other)
trajectories.assigned = read_excel("all_EUNIS_combinations_hand_assigned.xlsx", 
                                   sheet = 1) %>% 
  dplyr::select(Trajectory, Expert_System1, Expert_System2)

# As the file only contains habitat shifts from lower to higher alphabetical order
# the shifts from higher to lower order must be added, just with opposite 
# trajectories (i.e., succession becomes disturbance and vice versa)
trajectories.assigned2 = trajectories.assigned %>% 
  dplyr::select(Trajectory = 1, "Expert_System1" = 3, "Expert_System2" = 2) %>% 
  mutate(Trajectory = recode(Trajectory, 
                             "succession" = "disturbance", 
                             "disturbance" = "succession"))

trajectories.assigned = trajectories.assigned %>% 
  bind_rows(trajectories.assigned2)

rm(list = "trajectories.assigned2")

# load EUNIS habitat names (https://doi.org/10.1111/avsc.12519)
EUNIS.names = read_delim("EUNIS-ESy-2025-01-18_legendfromESYfile.txt",
                         delim = "\t") %>% 
  dplyr::select(code, name) %>% 
  # replace MA in EUNIS classification with M
  mutate(code = gsub("MA", "M", code))

# data overview
length(unique(dat.resurvey$`sPlot concept unified`)) # no. of taxa
nrow(dat.resurvey) # no. of species observations
length(unique(dat.resurvey$PlotObservationID)) # no. vegetation plots
length(unique(dat.resurvey$time_seriesID)) # no. of time series

# 2. calculate diversity indices ----------------------------------------------

# calculate compositional diversity
a = Sys.time()
resurvey.compositional.diversity = 
  dat.resurvey[,.(species_richness = calculate.species.richness(`sPlot concept unified`),
                  cover_sum = 1- exp(1)^(sum(log(1-(`Cover %` / 100)))),
                  shannon_diversity = calculate.shannon.diversity(`sPlot concept unified`, `Cover %`)),
               by = list(PlotObservationID, time_seriesID, `Date of recording`,
                         `ReSurvey project`, `ReSurvey site`, `ReSurvey plot`,
                         `Expert system`,  Longitude, Latitude)] # 4-6 min
Sys.time() - a

# calculate Pilou evenness 
resurvey.compositional.diversity$species_evenness = 
  resurvey.compositional.diversity$shannon_diversity / log(resurvey.compositional.diversity$species_richness)

# Calculate cover and no. of Red Dist species
# First, check if species are on European red list. If not, check if species are 
# on national Red List. Always use the more severe status.
resurvey.redlist = dat.resurvey %>% 
  distinct(PlotObservationID) %>% 
  left_join(dat.resurvey %>% 
              # join Red List information
              dplyr::select(PlotObservationID, Country, `sPlot concept unified`, `Cover %`) %>% 
              left_join(dat.redlist.europe %>%  
                          dplyr::select(Country, `sPlot concept unified`, Threat),
                        by = c("Country","sPlot concept unified")) %>%
              left_join(dat.redlist.national %>%  
                          dplyr::select(Country, `sPlot concept unified`, Threat_rl_NationalRL = Threat),
                        by = c("Country","sPlot concept unified")) %>% 
              # unify into one threat category (the highest threat level reported)
              mutate(Threat = ifelse(is.na(Threat), 0, Threat)) %>% 
              mutate(Threat_rl_NationalRL = ifelse(is.na(Threat_rl_NationalRL), 0, Threat_rl_NationalRL)) %>% 
              mutate(Threat = ifelse(Threat == 0 & Threat_rl_NationalRL == 1, 1, Threat)) %>% 
              # calculate number and cover of threatened species
              group_by(PlotObservationID, Threat, .drop = F) %>% 
              summarise(redlist_threatened_nr = length(`sPlot concept unified`),,
                        redlist_threatened_cover = 1- exp(1)^(sum(log(1-(`Cover %` / 100))))) %>% 
              filter(Threat == 1) %>% 
              dplyr::select(-Threat), 
            by = "PlotObservationID") %>% 
  # replace NAs in plots that do not have any threatened species
  mutate(redlist_threatened_nr = ifelse(is.na(redlist_threatened_nr), 0, redlist_threatened_nr),
         redlist_threatened_cover = ifelse(is.na(redlist_threatened_cover), 0, redlist_threatened_cover))

# calculate functional diversity
a = Sys.time()
resurvey.functional.diversity = dat.resurvey %>% 
  group_by(PlotObservationID) %>% 
  do(calculate.functional.diversity(.$`sPlot concept unified`, .$`Cover %`)) %>% 
  ungroup()
Sys.time() - a

# calculate phylogenetic diversity
a = Sys.time()
resurvey.phylogenetic.diversity = dat.resurvey %>% 
  mutate(`sPlot concept unified` = gsub(" ", "_", `sPlot concept unified`)) %>% 
  group_by(PlotObservationID) %>% 
  do(calculate.phylogenetic.diversity(.$`sPlot concept unified`, .$`Cover %`)) %>% 
  ungroup()
Sys.time() - a # ~50-60 min

# bind all diversity information together
resurvey.diversity = resurvey.compositional.diversity %>% 
  bind_cols(resurvey.functional.diversity %>% dplyr::select(- PlotObservationID),
            resurvey.phylogenetic.diversity %>%  dplyr::select(- PlotObservationID),
            resurvey.redlist %>%  dplyr::select(- PlotObservationID))

# save plot-level diversity values
fwrite(resurvey.diversity, 
       "plot_level_diversity.csv",
       sep = "\t")


# calculate decadal trends in species pool richness
# -> clunky, but does the trick as speed is not an issue
results.decadal.trends = list()
for(year.temp in seq(from = 1910, to = 2020, by = 10)){
  results.decadal.trends[[as.character(year.temp)]] = 
    get.decadal.changes.in.species.richness.with.hand.assigned.trajectories(
      dat.resurvey, year.temp, timespan.to.analyse = 10)}

# bind rows of nested lists
for(i in 1:length(results.decadal.trends)){
  for(j in 1:length(results.decadal.trends[[i]])){
    for(k in 1:length(results.decadal.trends[[i]][[j]])){
      results.decadal.trends[[i]][[j]][[k]] = bind_rows(results.decadal.trends[[i]][[j]][[k]])}
    results.decadal.trends[[i]][[j]] = bind_rows(results.decadal.trends[[i]][[j]])}
  results.decadal.trends[[i]] = bind_rows(results.decadal.trends[[i]])}

results.decadal.trends.based.on.hand.assigned.trajectories = bind_rows(results.decadal.trends)

# save decadal trends 
write_delim(results.decadal.trends.based.on.hand.assigned.trajectories, 
            "decadal_trends_in_pool_richness_EUNIS3.csv", 
            delim = "\t")

# 3. Prepare data for analyses --------------------------------------------
resurvey.diversity = fread("plot_level_diversity.csv", sep = "\t")

# sort, convert, correct data on plot-level diversity indices
resurvey.diversity = resurvey.diversity %>% 
  # sort from earliest to latest sampling dates
  arrange(`Date of recording`, time_seriesID) %>% 
  # convert Date format
  mutate(`Date of recording` = as.Date(`Date of recording`, tryFormats = c("%d.%m.%Y"))) %>% 
  mutate(year_of_recording = year(`Date of recording`)) %>% 
  mutate(`Date of recording` = as.character(`Date of recording`)) %>% 
  # correct wrong coordinate format
  mutate(corrected_lon_lat = ifelse((Longitude > 200 | Longitude < -200) & (Latitude > 200 | Latitude < -200), T, F)) %>%
  mutate(Latitude = ifelse(corrected_lon_lat == T, Latitude / 10000, Latitude),
         Longitude = ifelse(corrected_lon_lat == T, Longitude / 10000, Longitude)) %>%
  # replace MA in EUNIS classification with M
  mutate(`Expert system` = gsub("MA", "M", `Expert system`))  %>% 
  # determine common EUNIS levels for plots that were classified in multiple habitats
  rowwise() %>% 
  mutate(`Expert system` = ifelse(`Expert system` %in% c("", "~"), "~", `Expert system`)) %>% 
  mutate(`Expert system` = unify.EUNIS.levels(`Expert system`)) %>%
  ungroup() # 1 min

# For following functions, replace white spaces and minus signs in column names
names(resurvey.diversity) = gsub(" ", "_", names(resurvey.diversity))  
names(resurvey.diversity) = gsub("-", "_", names(resurvey.diversity))  

# define all diversity indices that should be analysed
all.metrics = c("cover_sum", "species_richness", "shannon_diversity", 
                "FDis", "FDiv", "FEve", "PD", "MPD", "MNTD",
                "redlist_threatened_nr", "redlist_threatened_cover")

# the same of the indices for annual percentage changes
all.rate.metrics = c("cover_sum", "species_richness", "shannon_diversity",
                     "FDis", "FDiv", "FEve", "PD", "MPD", "MNTD", 
                     "redlist_threatened_nr", "redlist_threatened_cover")

# provide a definition of full EUNIS level 1 habitat names
EUNIS1.rename = c("~" = "- unassigned", 
                  "M" = "MA - Marine habitats", 
                  "N" = "N - Coastal habitats",
                  "P" = "P - Inland waters",
                  "Q" = "Q - Mires and wetlands",
                  "R" = "R - Grasslands", 
                  "S" = "S - Shrublands and heathlands", 
                  "T" = "T - Forests", 
                  "U" = "U - Inland sparse vegetation", 
                  "V" = "V - Vegetated man-made habitats")

# define the colours for plotting
EUNIS1.colors = c("all" = "black", 
                  "- unassigned" = "#D3D3D3",
                  "MA - Marine habitats" = "#bf93d2",
                  "N - Coastal habitats" = "#9ecae1",
                  "P - Inland waters" = "#6083d0",
                  "Q - Mires and wetlands" = "#addd8e",
                  "R - Grasslands" = "#31a354", 
                  "S - Shrublands and heathlands" = "#ffffd4", 
                  "T - Forests" = "#fe9929", 
                  "U - Inland sparse vegetation" = "#eaa9a5", 
                  "V - Vegetated man-made habitats" = "#b7251b")

# 4a. calculate linear annual trends of absolute changes per time series -------

a = Sys.time()
for(metric.temp in all.metrics){
  
  # separately for each time_seriesID
  results.temp = resurvey.diversity %>% 
    dplyr::select(time_seriesID, div_value = all_of(metric.temp), Expert_system, year_of_recording) %>% 
    filter(!is.infinite(div_value)) %>% 
    drop_na() %>% 
    group_split(time_seriesID) %>% 
    lapply(FUN = analyse.and.compile.lms, metric.temp) %>% 
    bind_rows()
  
  write_delim(results.temp, 
              paste0("lm_absolute_results/lm_", metric.temp, ".csv"),
              delim = "\t")
  print(paste0("done: ", metric.temp, " - ",  round(difftime(Sys.time(), a, units='mins'), 2), " min"))
}

rm("results.temp"); gc()

# 4b. calculate annual percentage changes per time series and ------------------
# relative change is calculated by using lm(log(x) ~ year_of_recording)
# followed by a back-transformation of the coefficient estimate as 
# = exp(coef) - 1, which gives the annual percentage change

a = Sys.time()
for(metric.temp in all.rate.metrics){
  
  # separate by time series id
  results.temp = resurvey.diversity %>% 
    dplyr::select(time_seriesID, div_value = all_of(metric.temp), Expert_system, year_of_recording) %>% 
    filter(!is.infinite(div_value))
  
  # redlist metrics must be added a constant of 0.5 since these can be zero
  # which is half of the minimum value (0.001 would lead to wrong coefficient estimates)
  if(metric.temp %in% c("redlist_threatened_nr", "redlist_threatened_cover")){
    results.temp$div_value = results.temp$div_value + 0.5}
  
  results.temp = results.temp %>% 
    drop_na() %>% 
    group_split(time_seriesID) %>% 
    lapply(FUN = analyse.and.compile.lms.of.rate.change, metric.temp) %>% 
    bind_rows()
  
  write_delim(results.temp, 
              paste0("lm_rate_results/lm_", metric.temp, ".csv"),
              delim = "\t")
  print(paste0("done: ", metric.temp, " - ",  round(difftime(Sys.time(), a, units='mins'), 2), " min"))
}

rm("results.temp"); gc()

# 5. re-organize results of linear and rate trends analyses --------------------
lm.model.results = list()
lm.rate.model.results = list()

for(i in 1:length(all.metrics)){
  metric.temp = all.metrics[i]
  lm.model.results[[metric.temp]] = read_delim(paste0("lm_absolute_results/lm_",metric.temp, ".csv"),
                                               delim = "\t")
  print(paste0(metric.temp, " - ", nrow(lm.model.results[[metric.temp]]), " rows"))}

for(i in 1:length(all.rate.metrics)){
  metric.temp = all.rate.metrics[i]
  lm.rate.model.results[[metric.temp]] = read_delim(paste0("lm_rate_results/lm_",metric.temp, ".csv"),
                                                    delim = "\t")
  print(paste0(metric.temp, " - ", nrow(lm.rate.model.results[[metric.temp]]), " rows"))}

# extract slopes  
all.lm.tidy.fits = lm.model.results %>% 
  bind_rows() %>% 
  filter(n > 1) %>% 
  mutate(weight_value = log(n),
         time_span = year_end - year_start)

all.rate.lm.tidy.fits = lm.rate.model.results %>% 
  bind_rows() %>% 
  filter(n > 1) %>% 
  mutate(weight_value = log(n),
         time_span = year_end - year_start)

# free memory space
remove("lm.model.results")
remove("lm.rate.model.results")
gc()

# remove time series that have only data within one year
all.lm.tidy.fits = all.lm.tidy.fits %>% 
  filter(year_end != year_start)

all.rate.lm.tidy.fits = all.rate.lm.tidy.fits %>% 
  filter(year_end != year_start)

# correct auxiliary EUNIS classes Qa + Qb -> Q | Sa + Sb -> S
# which are auxilliary classes
all.lm.tidy.fits = all.lm.tidy.fits %>% 
  mutate(EUNIS_start = ifelse(EUNIS_start %in% c("Qa", "Qb"), "Q", EUNIS_start),
         EUNIS_start = ifelse(EUNIS_start %in% c("Sa", "Sb"), "S", EUNIS_start),
         EUNIS_end = ifelse(EUNIS_end %in% c("Qa", "Qb"), "Q", EUNIS_end),
         EUNIS_end = ifelse(EUNIS_end %in% c("Sa", "Sb"), "S", EUNIS_end))

all.rate.lm.tidy.fits = all.rate.lm.tidy.fits %>% 
  mutate(EUNIS_start = ifelse(EUNIS_start %in% c("Qa", "Qb"), "Q", EUNIS_start),
         EUNIS_start = ifelse(EUNIS_start %in% c("Sa", "Sb"), "S", EUNIS_start),
         EUNIS_end = ifelse(EUNIS_end %in% c("Qa", "Qb"), "Q", EUNIS_end),
         EUNIS_end = ifelse(EUNIS_end %in% c("Sa", "Sb"), "S", EUNIS_end))

# separate the levels of the EUNIS classification
all.lm.tidy.fits = all.lm.tidy.fits %>% 
  mutate(EUNIS1_start = substr(EUNIS_start, start = 1, stop = 1),
         EUNIS2_start = substr(EUNIS_start, start = 1, stop = 2),
         EUNIS3_start = substr(EUNIS_start, start = 1, stop = 3)) %>% 
  mutate(EUNIS1_end = substr(EUNIS_end, start = 1, stop = 1),
         EUNIS2_end = substr(EUNIS_end, start = 1, stop = 2),
         EUNIS3_end = substr(EUNIS_end, start = 1, stop = 3)) %>% 
  # remove classifications that are not on the according level
  mutate(EUNIS2_start = ifelse(nchar(EUNIS2_start) < 2, NA, EUNIS2_start),
         EUNIS3_start = ifelse(nchar(EUNIS3_start) < 3, NA, EUNIS3_start)) %>% 
  mutate(EUNIS2_end = ifelse(nchar(EUNIS2_end) < 2, NA, EUNIS2_end),
         EUNIS3_end = ifelse(nchar(EUNIS3_end) < 3, NA, EUNIS3_end))

all.rate.lm.tidy.fits = all.rate.lm.tidy.fits %>% 
  mutate(EUNIS1_start = substr(EUNIS_start, start = 1, stop = 1),
         EUNIS2_start = substr(EUNIS_start, start = 1, stop = 2),
         EUNIS3_start = substr(EUNIS_start, start = 1, stop = 3)) %>% 
  mutate(EUNIS1_end = substr(EUNIS_end, start = 1, stop = 1),
         EUNIS2_end = substr(EUNIS_end, start = 1, stop = 2),
         EUNIS3_end = substr(EUNIS_end, start = 1, stop = 3)) %>% 
  # remove classifications that are not on the according level
  mutate(EUNIS2_start = ifelse(nchar(EUNIS2_start) < 2, NA, EUNIS2_start),
         EUNIS3_start = ifelse(nchar(EUNIS3_start) < 3, NA, EUNIS3_start)) %>% 
  mutate(EUNIS2_end = ifelse(nchar(EUNIS2_end) < 2, NA, EUNIS2_end),
         EUNIS3_end = ifelse(nchar(EUNIS3_end) < 3, NA, EUNIS3_end))

# rename level 1 habitat types
for(i in 1:length(EUNIS1.rename)){
  all.lm.tidy.fits$EUNIS1_end[all.lm.tidy.fits$EUNIS1_end == names(EUNIS1.rename)[i]] = EUNIS1.rename[i]
  all.lm.tidy.fits$EUNIS1_start[all.lm.tidy.fits$EUNIS1_start == names(EUNIS1.rename)[i]] = EUNIS1.rename[i]}

for(i in 1:length(EUNIS1.rename)){
  all.rate.lm.tidy.fits$EUNIS1_start[all.rate.lm.tidy.fits$EUNIS1_start == names(EUNIS1.rename)[i]] = EUNIS1.rename[i]
  all.rate.lm.tidy.fits$EUNIS1_end[all.rate.lm.tidy.fits$EUNIS1_end == names(EUNIS1.rename)[i]] = EUNIS1.rename[i]}

# rename unclassified EUNIS level 2 and level 3 habitats
all.lm.tidy.fits = all.lm.tidy.fits %>% 
  mutate(EUNIS2_start = ifelse(is.na(EUNIS2_start) | EUNIS2_start %in% c("", "~"), 
                               "- unassigned", EUNIS2_start),
         EUNIS2_end = ifelse(is.na(EUNIS2_end) | EUNIS2_end %in% c("", "~"), 
                             "- unassigned", EUNIS2_end),
         EUNIS3_start = ifelse(is.na(EUNIS3_start) | EUNIS3_start %in% c("", "~"), 
                               "- unassigned", EUNIS3_start),
         EUNIS3_end = ifelse(is.na(EUNIS3_end) | EUNIS3_end %in% c("", "~"), 
                             "- unassigned", EUNIS3_end))

all.rate.lm.tidy.fits = all.rate.lm.tidy.fits %>% 
  mutate(EUNIS2_start = ifelse(is.na(EUNIS2_start) | EUNIS2_start %in% c("", "~"), 
                               "- unassigned", EUNIS2_start),
         EUNIS2_end = ifelse(is.na(EUNIS2_end) | EUNIS2_end %in% c("", "~"), 
                             "- unassigned", EUNIS2_end),
         EUNIS3_start = ifelse(is.na(EUNIS3_start) | EUNIS3_start %in% c("", "~"), 
                               "- unassigned", EUNIS3_start),
         EUNIS3_end = ifelse(is.na(EUNIS3_end) | EUNIS3_end %in% c("", "~"), 
                             "- unassigned", EUNIS3_end))

# attach information on EUNIS habitat trajectories
all.lm.tidy.fits = all.lm.tidy.fits %>% 
  # match with hand-assigned habitat-change trajectories
  left_join(EUNIS.trajectories,
            by = join_by("EUNIS_start" == "Expert_System1", "EUNIS_end" == "Expert_System2")) %>% 
  
  # when EUNIS3 is the same at start and end -> "stable"
  mutate(Trajectory = ifelse(is.na(Trajectory) & EUNIS3_start != "- unassigned" & EUNIS3_start == EUNIS3_end, "stable", Trajectory)) %>% 
  
  # when there is an "~" in EUNIS level 1 start AND end -> "unassigned"
  mutate(Trajectory = ifelse(is.na(Trajectory) & EUNIS_start == "~" & EUNIS_end == "~", "unassigned", Trajectory)) %>% 
  
  # when there is an "~" in either EUNIS level 1 start OR end -> "other"
  mutate(Trajectory = ifelse(is.na(Trajectory) & (EUNIS_start == "~" | EUNIS_end == "~"), "other", Trajectory)) %>% 
  
  # when only EUNIS level 2 for start and end available -> "unassigned
  mutate(Trajectory = ifelse(is.na(Trajectory) & (EUNIS2_start == EUNIS2_end), "unassigned", Trajectory)) %>% 
  
  # when only EUNIS level 1 for start and end available -> "unassigned"
  mutate(Trajectory = ifelse(is.na(Trajectory) & (EUNIS1_start == EUNIS1_end), "unassigned", Trajectory))

all.rate.lm.tidy.fits = all.rate.lm.tidy.fits %>% 
  # match with hand-assigned habitat-change trajectories
  left_join(EUNIS.trajectories,
            by = join_by("EUNIS_start" == "Expert_System1", "EUNIS_end" == "Expert_System2")) %>% 
  
  # when EUNIS3 is the same at start and end -> "stable"
  mutate(Trajectory = ifelse(is.na(Trajectory) & EUNIS3_start != "- unassigned" & EUNIS3_start == EUNIS3_end, "stable", Trajectory)) %>% 
  
  # when there is an "~" in EUNIS level 1 start AND end -> "unassigned"
  mutate(Trajectory = ifelse(is.na(Trajectory) & EUNIS_start == "~" & EUNIS_end == "~", "unassigned", Trajectory)) %>% 
  
  # when there is an "~" in either EUNIS level 1 start OR end -> "other"
  mutate(Trajectory = ifelse(is.na(Trajectory) & (EUNIS_start == "~" | EUNIS_end == "~"), "other", Trajectory)) %>% 
  
  # when only EUNIS level 2 for start and end available -> "unassigned
  mutate(Trajectory = ifelse(is.na(Trajectory) & (EUNIS2_start == EUNIS2_end), "unassigned", Trajectory)) %>% 
  
  # when only EUNIS level 1 for start and end available -> "unassigned"
  mutate(Trajectory = ifelse(is.na(Trajectory) & (EUNIS1_start == EUNIS1_end), "unassigned", Trajectory))

# save data
saveRDS(all.lm.tidy.fits, 
        "all_lm_tidy_fits.RDS")
saveRDS(all.rate.lm.tidy.fits, 
        "all_rate_lm_tidy_fits.RDS")
gc()

# H1: Is there a general decline in local plant diversity (absolute values) ----
# to analyse annual percentage changes, load file "all_rate_lm_tidy_fits.RDS"
all.lm.tidy.fits = readRDS("all_lm_tidy_fits.RDS") %>% 
  # rename diversity indices
  mutate(metric = recode(metric,
                         "cover_sum" = "Cover", "redlist_threatened_cover" = "Cover of\nred list species", "redlist_threatened_nr" = "Number of\nred list species",
                         "species_richness" = "Species\nrichness", "shannon_diversity" = "Shannon\ndiversity",
                         "FDis" = "Functional\ndispersion", "FDiv" = "Functional\ndivergence", "FEve" = "Functional\nevenness", 
                         "PD" = "Faith\nphylogenetic diversity", "MPD" = "Mean pairwise\nphylogenetic distance", "MNTD" = "Mean nearest\ntaxon distance")) 

# calculate grand mean trends across all time series
# using linear mixed-effect model with random effect for EUNIS level 1 habitat type
grand.means.all = all.lm.tidy.fits %>% 
  drop_na(slope) %>% 
  dplyr::select(slope, weight_value, metric, EUNIS1_start) %>% 
  nest(lme_data = c(slope, weight_value, EUNIS1_start)) %>% 
  mutate(lme_fit = purrr::map(lme_data, ~ lme4::lmer(slope ~ 1 + (1|EUNIS1_start), 
                                                     weight = weight_value, data = .x))) %>% 
  mutate(tidy_fit = purrr::map(lme_fit, broom.mixed::tidy, conf.int = T, conf.method = "Wald")) 

# collect results in tibble
grand.means.all = grand.means.all$tidy_fit %>%
  bind_rows() %>% 
  filter(effect == "fixed") %>% 
  mutate(metric = grand.means.all$metric) %>% 
  dplyr::select(metric, estimate, conf.low,conf.high) %>% 
  mutate(estimate_round = round(estimate, 3))

# H2: variance partitioning ----------------------------------------------------
# to analyse annual percentage changes, load file "all_rate_lm_tidy_fits.RDS"
all.lm.tidy.fits = readRDS("all_lm_tidy_fits.RDS") %>% 
  mutate(metric = recode(metric,
                         "cover_sum" = "Cover", "redlist_threatened_cover" = "Cover of\nred list species", "redlist_threatened_nr" = "Number of\nred list species",
                         "species_richness" = "Species\nrichness", "shannon_diversity" = "Shannon\ndiversity",
                         "FDis" = "Functional\ndispersion", "FDiv" = "Functional\ndivergence", "FEve" = "Functional\nevenness", 
                         "PD" = "Faith\nphylogenetic diversity", "MPD" = "Mean pairwise\nphylogenetic distance", "MNTD" = "Mean nearest\ntaxon distance")) 

# filter data
data.for.varpart = all.lm.tidy.fits %>% 
  # EUNIS3 must be available
  filter(EUNIS3_start != "- unassigned") %>% 
  # Trajectory must be clear
  filter(Trajectory %in% c("stable", "succession", "disturbance")) %>% 
  # at least 10 time series per trajectory type
  group_by(metric, EUNIS3_start, Trajectory) %>% 
  mutate(trajectory_n = length(time_seriesID)) %>% 
  ungroup() %>% 
  filter(trajectory_n >= 10) %>%
  # at least 2 different trajectories per EUNIS3 habitat type
  group_by(metric, EUNIS3_start) %>% 
  mutate(n_trajectories = length(unique(Trajectory))) %>% 
  ungroup() %>% 
  filter(metric >= 2) %>% 
  # at least two EUNIS3_start per EUNIS1_start
  group_by(metric, EUNIS1_start) %>% 
  mutate(EUNIS3_types = length(unique(EUNIS3_start))) %>% 
  ungroup() %>% 
  filter(EUNIS3_types >= 2)

# conduct variance partitioning within each EUNIS1 habitat type
results.temp = list()

for(metric.temp in unique(data.for.varpart$metric)){
  
  results.temp[[metric.temp]] = list()
  
  dat.temp = data.for.varpart %>%  
    filter(metric == metric.temp)  
  
  for(EUNIS1.temp in unique(dat.temp$EUNIS1_start)){
    
    dat.temp.temp = dat.temp %>% 
      filter(EUNIS1_start == EUNIS1.temp)
    
    # use hierarchical partitioning from Chevan and Sutherland (1991)
    library(ghp)
    full_model = lm(slope ~ EUNIS3_start + Trajectory + year_end, data = dat.temp.temp)
    part.temp = ghp(depname = "slope", dat.temp.temp[c("slope", "EUNIS3_start", "Trajectory", "year_end")], 
                    method = "lm", gof = "r.squared")
    
    results.temp[[metric.temp]][[EUNIS1.temp]] = part.temp$results %>% 
      dplyr::select(var, indep_effects) %>% 
      pivot_wider(names_from = "var", 
                  values_from = "indep_effects") %>% 
      mutate(Unexplained = 1 - summary(full_model)$r.squared) %>% 
      mutate(EUNIS1_start = EUNIS1.temp, metric = metric.temp,
             p_value = summary(full_model)$fstatistic %>% {unname(pf(.[1],.[2],.[3],lower.tail=F))}) %>% 
      mutate(model_sign = ifelse(p_value < 0.05, "significant", "non-significant"), 
             sign_sign = ifelse(p_value < 0.05, "*", " "))
  }
  print(paste0("done: ", metric.temp))
}


# flatten the obtained nested list
for(i in 1:length(results.temp)){
  results.temp[[i]] = bind_rows(results.temp[[i]])}
results.temp = bind_rows(results.temp)

# H2: Specificity of diversity trends (absolute values) ------------------------
# to analyse annual percentage changes, load file "all_rate_lm_tidy_fits.RDS"
all.lm.tidy.fits = readRDS("all_lm_tidy_fits.RDS") %>% 
  mutate(metric = recode(metric,
                         "cover_sum" = "Cover", "redlist_threatened_cover" = "Cover of\nred list species", "redlist_threatened_nr" = "Number of\nred list species",
                         "species_richness" = "Species\nrichness", "shannon_diversity" = "Shannon\ndiversity",
                         "FDis" = "Functional\ndispersion", "FDiv" = "Functional\ndivergence", "FEve" = "Functional\nevenness", 
                         "PD" = "Faith\nphylogenetic diversity", "MPD" = "Mean pairwise\nphylogenetic distance", "MNTD" = "Mean nearest\ntaxon distance")) 

# Calculate group-level mean trends within EUNIS1 habitats
EUNIS1.grand.means = all.lm.tidy.fits %>% 
  drop_na(slope) %>% 
  filter(EUNIS1_start != "- unassigned") %>% 
  dplyr::select(slope, weight_value, metric, EUNIS1_start) %>% 
  nest(lm_data = c(slope, weight_value)) %>% 
  mutate(lm_fit = map(lm_data, ~ lm(slope ~ 1, 
                                    weight = weight_value, data = .x)),
         tidy_fit = map(lm_fit, broom.mixed::tidy, conf.int = T, conf.method = "Wald"),
         n = map_dbl(lm_data, nrow)) # 1 min

# Calculate group-level mean trends within EUNIS1 habitats and habitat-change trajectories
EUNIS1.means.per.trajectory = all.lm.tidy.fits %>% 
  drop_na(slope) %>% 
  filter(EUNIS1_start != "- unassigned") %>% 
  dplyr::select(slope, weight_value, metric, EUNIS1_start, Trajectory) %>% 
  nest(lm_data = c(slope, weight_value)) %>% 
  mutate(lm_fit = map(lm_data, ~ lm(slope ~ 1, 
                                    weight = weight_value, data = .x)),
         tidy_fit = map(lm_fit, broom.mixed::tidy, conf.int = T, conf.method = "Wald"),
         n = map_dbl(lm_data, nrow)) # 1 min

# collect results
EUNIS1.grand.means = EUNIS1.grand.means$tidy_fit %>%
  bind_rows() %>% 
  mutate(metric = EUNIS1.grand.means$metric,
         EUNIS1_start = EUNIS1.grand.means$EUNIS1_start,
         n = EUNIS1.grand.means$n)

EUNIS1.means.per.trajectory = EUNIS1.means.per.trajectory$tidy_fit %>%
  bind_rows() %>% 
  mutate(metric = EUNIS1.means.per.trajectory$metric,
         EUNIS1_start = EUNIS1.means.per.trajectory$EUNIS1_start, 
         Trajectory = EUNIS1.means.per.trajectory$Trajectory,
         n = EUNIS1.means.per.trajectory$n)


# H3: Decadal trends in species pool richness ----------------------------------
dat.decadal.richness = read_delim("decadal_trends_in_pool_richness_EUNIS3.csv", 
                                  delim = "\t")

# only include decadal trends with >= 10 time series available
dat.decadal.richness = dat.decadal.richness %>% 
  filter(n >= 10) %>% 
  group_by(trajectory, EUNIS1) %>% 
  mutate(nr_of_decadal_comparisons = length(decade1)) %>% 
  ungroup() %>% 
  filter(nr_of_decadal_comparisons >= 2) %>% 
  # make categories for sample size
  mutate(n_cut = ifelse(n <= 20, "≤ 20", 
                        ifelse(n > 200, "> 200", "≤ 200"))) %>% 
  mutate(n_cut = factor(n_cut, levels = c("≤ 20", "≤ 200", "> 200")))

# calculate decadal linear trends
decadal.trends = dat.decadal.richness %>% 
  dplyr::select(-local_richness1, - local_richness2) %>% 
  mutate(decade1 = decade1 + 5, 
         decade2 = decade2 + 5,
         pool_delta = pool_richness2 - pool_richness1,
         pool_delta_per_decade = (pool_richness2 - pool_richness1) / ((decade2 - decade1)/10), 
         pool_perc_change = (pool_richness2 - pool_richness1)/ pool_richness1 * 100,
         timespan = decade2 - decade1) %>% 
  mutate(pool_trend = ifelse(pool_delta > 0, "positive", 
                             ifelse(pool_delta < 0, "negative", "0"))) %>% 
  mutate(pool_trend = factor(pool_trend, levels = c("positive", "0", "negative"))) %>% 
  mutate(pool_trend_numeric = 0) %>% 
  mutate(pool_trend_numeric = ifelse(pool_trend == "positive", 1, -1))

# for pool richness: calculate t- and binomial tests of trends
list.of.decadal.tests = list()

# conduct t-test and binomial test for each decadal comparison
for(trajectory.temp in unique(decadal.trends$trajectory)){
  
  list.of.decadal.tests[[trajectory.temp]] = list()
  
  dat.temp = decadal.trends %>% 
    filter(trajectory == trajectory.temp)
  
  # for each habitat
  for(eunis.temp in unique(dat.temp$EUNIS1)){
    
    dat.temp.temp = filter(dat.temp, EUNIS1 == eunis.temp)
    
    # tests for pool richness
    test.pool.t.temp = t.test(dat.temp.temp$pool_delta_per_decade, mu = 0)
    test.pool.binomial.temp = 
      binom.test(x = dat.temp.temp %>% 
                   filter(pool_trend == "positive") %>% 
                   nrow(), 
                 n = dat.temp.temp %>% 
                   filter(pool_trend %in% c("positive", "negative")) %>% 
                   nrow())
    
    # collect results
    list.of.decadal.tests[[trajectory.temp]][[eunis.temp]] = tibble(
      "trajectory" = trajectory.temp,
      "EUNIS1" = eunis.temp,
      "t_test_mean" = test.pool.t.temp$estimate,
      "t_test_p" = test.pool.t.temp$p.value,
      "binom_test_proportion_positive" = test.pool.binomial.temp$estimate,
      "binom_test_p" = test.pool.binomial.temp$p.value,
      xpos = -Inf, ypos = Inf, hjustvar = -0.02, vjustvar = 1.2)
  }
}

# flatten test results
for(i in 1:length(list.of.decadal.tests)){
  list.of.decadal.tests[[i]] = bind_rows(list.of.decadal.tests[[i]])}

# flatten and assign significance labels
list.of.decadal.tests = list.of.decadal.tests %>% 
  bind_rows() %>% 
  mutate(t_test_sign = ifelse(t_test_p < 0.05, "*", ""),
         binom_test_sign = ifelse(binom_test_p < 0.05, "*", ""))
