# find common habitat for plots with multiple EUNIS habitats assigned
unify.EUNIS.levels = function(eunis.temp){
  eunis.sep = str_split(eunis.temp, pattern = ",")[[1]]
  eunis.sep = gsub("!", "", eunis.sep)
  if(length(unique(substr(eunis.sep, 1, 3))) == 1){
    return(substr(eunis.sep, 1,3)[1])}
  if(length(unique(substr(eunis.sep, 1, 2))) == 1){
    return(substr(eunis.sep, 1,2)[1])}
  if(length(unique(substr(eunis.sep, 1, 1))) == 1){
    return(substr(eunis.sep, 1,1)[1])}
  return("~")
}

# calculate species richness -------
calculate.species.richness = function(taxon.temp){
  return(length(unique(taxon.temp)))}


# calculate Shannon diversity ---------
calculate.shannon.diversity = function(taxon.temp, cover.temp){
  
  if(sum(cover.temp,na.rm = T) > 0){
    return(data.table(taxon = taxon.temp, cover = cover.temp) %>% 
             reshape2::dcast(1 ~ taxon,
                             value.var = "cover",
                             fun.aggregate = sum) %>% 
             diversity(index = "shannon"))
  }else{
    return(as.double(NA))}
}

# calculate functional diversity indices -----------
calculate.functional.diversity = function(taxon.temp, cover.temp){
  
  # order dat.temp by taxon name
  dat.temp = data.table(taxon = taxon.temp,
                        cover_perc = cover.temp)[order(taxon)]
  
  # sum up cover values from different layers
  dat.temp = dat.temp[,.(cover_perc = sum(cover_perc, na.rm = T)), by = taxon]
  
  # for presence absence data, replace zero cover values with a value of 1
  if(sum(dat.temp$cover_perc, na.rm = T) == 0){
    dat.temp$cover_perc = 1}
  
  # join trait data
  traits.temp = dat.try.scaled[dat.temp, on = .(Taxon_name = taxon), nomatch = 0]
  
  # remove species that have no trait data
  traits.temp.no.na = na.omit(traits.temp)
  
  # do not calculate any div index if trait data available for less than 80% of the summed cover or 
  if(sum(traits.temp.no.na$cover_perc) / sum(dat.temp$cover_perc, na.rm = T) < 0.8){
    return(data.frame("FDis" = as.double(NA), "FDiv" = as.double(NA), "FEve" = as.double(NA), 
             "FRic" = as.double(NA), "qual.FRic" = as.double(NA), "RaoQ" = as.double(NA)))}
  
  # calculate only CWM/CWV when only one species with trait data available
  if(nrow(traits.temp.no.na) <= 2){
    return(data.frame("FDis" = as.double(NA), "FDiv" = as.double(NA), "FEve" = as.double(NA), 
             "FRic" = as.double(NA), "qual.FRic" = as.double(NA), "RaoQ" = as.double(NA)))}
  
  # pepare species and trait data
  x.temp = data.frame(traits.temp.no.na[, - "cover_perc"]) %>% 
    column_to_rownames("Taxon_name")
  
  a.temp = data.table::dcast(traits.temp.no.na[,c("Taxon_name", "cover_perc")],
                             1 ~ Taxon_name,
                             value.var = "cover_perc")[,- 1] 
  rownames(a.temp) = "site"
  
  # calculate functional diversity indices
  fdist.temp = dbFD(x = x.temp,
                    a = a.temp,
                    stand.x = F,
                    calc.CWM = F,
                    messages = F,
                    m = 2)
  
  # return results in vector form (for speed purpose)
  return(data.frame("FDis" = ifelse(length(fdist.temp$FDis) == 0, as.double(NA), fdist.temp$FDis),
           "FDiv" = ifelse(length(fdist.temp$FDiv) == 0, as.double(NA), fdist.temp$FDiv),
           "FEve" = ifelse(length(fdist.temp$FEve) == 0, as.double(NA), fdist.temp$FEve),
           "FRic" = ifelse(length(fdist.temp$FRic) == 0, as.double(NA), fdist.temp$FRic),
           "qual.FRic" = ifelse(length(fdist.temp$qual.FRic) == 0, as.double(NA), fdist.temp$qual.FRic),
           "RaoQ" = ifelse(length(fdist.temp$RaoQ) == 0, as.double(NA), fdist.temp$RaoQ)))
}

######################################################
calculate.phylogenetic.diversity = function(taxon.temp, cover.temp){
  
  if(all(taxon.temp %in% phylotree$tip.label)){
    
    # arrange data to wide format matrix
    dat.temp = matrix(cover.temp, nrow = 1)
    colnames(dat.temp) = taxon.temp
    
    # prune tree
    tree.temp = phylo.prune(taxon.temp, phylotree)
    
    # calculate phylogenetic indices
    return(data.frame("PD" = picante::pd(dat.temp, tree.temp, include.root = F)$PD,
                      "MPD" = picante::mpd(dat.temp, cophenetic(tree.temp), abundance.weighted = T),
                      "MNTD" = picante::mntd(dat.temp, cophenetic(tree.temp), abundance.weighted = T)))
    
    # or return NAs
  }else{return(data.frame("PD" =  NA, "MPD" = NA, "MNTD" = NA))}
}

#######################################################
get.decadal.changes.in.species.richness.with.hand.assigned.trajectories = 
  function(dat.resurvey, year.temp, timespan.to.analyse){
    
    results.temp = list()
    
    for(timespan.temp in timespan.to.analyse * 1:20){
      
      results.temp[[as.character(timespan.temp)]] = list("all" = list(),
                                                         "EUNIS2_stable" = list(),
                                                         "EUNIS2_unstable" = list())
      
      # filter for decades to analyse
      dat.resurvey1 = dat.resurvey %>% 
        filter(`Year of recording` %in% (year.temp + (0:(timespan.to.analyse-1))))
      
      dat.resurvey2 = dat.resurvey %>% 
        filter(`Year of recording` %in% ((year.temp + timespan.temp) + (0:(timespan.to.analyse-1))))
      
      # filter only plots with assigned EUNIS 3 class
      dat.resurvey1 = dat.resurvey1 %>% 
        mutate(EUNIS1 = substr(`Expert system`, 1, 1),
               EUNIS3 = substr(`Expert system`, 1, 3),) %>% 
        mutate(EUNIS3 = ifelse(is.na(EUNIS3), "~", EUNIS3))
      
      dat.resurvey2 = dat.resurvey2 %>% 
        mutate(EUNIS1 = substr(`Expert system`, 1, 1),
               EUNIS3 = substr(`Expert system`, 1, 3)) %>% 
        mutate(EUNIS3 = ifelse(is.na(EUNIS3), "~", EUNIS3)) %>% 
        filter(nchar(EUNIS3) == 3)
      
      # filter only time series with single observations 
      # i.e. one PlotObservationID ReSurvey project x site x plot x date
      # omit N-to-N relationships
      dat.resurvey1 = dat.resurvey1 %>% 
        group_by(`ReSurvey project`, `ReSurvey site`, `ReSurvey plot`, `Date of recording`) %>% 
        mutate(n_PlotObservationID = length(unique(PlotObservationID))) %>% 
        ungroup() %>% 
        filter(n_PlotObservationID == 1)
      
      dat.resurvey2 = dat.resurvey2 %>% 
        group_by(`ReSurvey project`, `ReSurvey site`, `ReSurvey plot`, `Date of recording`) %>% 
        mutate(n_PlotObservationID = length(unique(PlotObservationID))) %>% 
        ungroup() %>% 
        filter(n_PlotObservationID == 1)
      
      # filter only times series that were sampled in decade1 and decade2
      dat.resurvey1 = dat.resurvey1 %>% 
        filter(time_seriesID %in% dat.resurvey2$time_seriesID)
      
      dat.resurvey2 = dat.resurvey2 %>% 
        filter(time_seriesID %in% dat.resurvey1$time_seriesID)
      
     # check if there is anything to analyse
      if(nrow(dat.resurvey1) > 0){
        
        # only first sampling in decade1 and last sampling in decade2
        dat.resurvey1 = dat.resurvey1 %>% 
          group_by(time_seriesID) %>% 
          mutate(earliest_Date = min(`Date of recording`)) %>% 
          filter(`Date of recording` == earliest_Date) %>% 
          ungroup()
        
        dat.resurvey2 = dat.resurvey2 %>% 
          group_by(time_seriesID) %>% 
          mutate(latest_Date = max(`Date of recording`)) %>% 
          filter(`Date of recording` == latest_Date) %>% 
          ungroup()
        
        # assign habitat shift trajectories
        # stable, succession, disturbance, other, or unassigned
        habitat.shift.trajectory.classification = dat.resurvey1 %>% 
          dplyr::select(time_seriesID, EUNIS_start = EUNIS3) %>% 
          distinct() %>% 
          full_join(dat.resurvey2 %>%  
                      dplyr::select(time_seriesID, EUNIS_end = EUNIS3) %>% 
                      distinct(),
                    by = "time_seriesID") %>% 
          left_join(trajectories.assigned %>% 
                      dplyr::select(EUNIS_trajectory = Trajectory, EUNIS_start = Expert_System1, EUNIS_end = Expert_System2),
                    by = c("EUNIS_start" = "EUNIS_start", "EUNIS_end" = "EUNIS_end"))
        
        # those trajectories that are not assigned by hand can either be
        # stable = same level 3 habitat type
        # undetermined = unassigned habitat type or same level 1 or level 2 habitat type
        habitat.shift.trajectory.classification = habitat.shift.trajectory.classification %>% 
          mutate(nr_digits_start =  nchar(gsub("[^0-9]+", "", EUNIS_start)),
                 nr_digits_end =  nchar(gsub("[^0-9]+", "", EUNIS_end))) %>% 
          mutate(EUNIS3_for_start_and_end = ifelse(nr_digits_start == 2 & nr_digits_end == 2,
                                                  "yes", "no")) %>% 
          mutate(EUNIS_trajectory = ifelse(is.na(EUNIS_trajectory) & 
                                             EUNIS3_for_start_and_end == "yes" & 
                                             EUNIS_start == EUNIS_end, "stable", EUNIS_trajectory)) %>% 
          mutate(EUNIS_trajectory = ifelse(is.na(EUNIS_trajectory) & 
                                             EUNIS3_for_start_and_end == "no", "unassigned", EUNIS_trajectory)) %>% 
          dplyr::select(- nr_digits_start , - nr_digits_end, - EUNIS3_for_start_and_end)
        
        # give output if a certain trajectory has not yet been classified
        if(any(is.na(habitat.shift.trajectory.classification$EUNIS_trajectory))){
          print("missing trajectory assignment")
          
          to.check.temp = habitat.shift.trajectory.classification %>% 
            filter(is.na(EUNIS_trajectory))  %>% 
            dplyr::select(EUNIS_start, EUNIS_end) %>% 
            distinct(EUNIS_start, EUNIS_end)
            
          to.check = to.check[[length(to.check) + 1]] = to.check.temp
          assign("to.check", to.check, envir = globalenv())
          just.throw.and.error.already # create error
        }
        
        time.seriesID.stable = habitat.shift.trajectory.classification %>% 
          filter(EUNIS_trajectory == "stable") %>%  
          pull(time_seriesID)
        
        time.seriesID.succession = habitat.shift.trajectory.classification %>% 
          filter(EUNIS_trajectory == "succession") %>%  
          pull(time_seriesID)
        
        time.seriesID.disturbance = habitat.shift.trajectory.classification %>% 
          filter(EUNIS_trajectory == "disturbance") %>%  
          pull(time_seriesID)
        
        time.seriesID.other = habitat.shift.trajectory.classification %>% 
          filter(EUNIS_trajectory == "other") %>%  
          pull(time_seriesID)
        
        time.seriesID.unassigned = habitat.shift.trajectory.classification %>% 
          filter(EUNIS_trajectory == "unassigned") %>%  
          pull(time_seriesID)
        
        dat.resurvey1.stable = as.data.table(dat.resurvey1)[time_seriesID %in% time.seriesID.stable]
        dat.resurvey2.stable = as.data.table(dat.resurvey2)[time_seriesID %in% time.seriesID.stable]
        
        dat.resurvey1.succession = as.data.table(dat.resurvey1)[time_seriesID %in% time.seriesID.succession]
        dat.resurvey2.succession = as.data.table(dat.resurvey2)[time_seriesID %in% time.seriesID.succession]
        
        dat.resurvey1.disturbance = as.data.table(dat.resurvey1)[time_seriesID %in% time.seriesID.disturbance]
        dat.resurvey2.disturbance = as.data.table(dat.resurvey2)[time_seriesID %in% time.seriesID.disturbance]
        
        dat.resurvey1.other = as.data.table(dat.resurvey1)[time_seriesID %in% time.seriesID.other]
        dat.resurvey2.other = as.data.table(dat.resurvey2)[time_seriesID %in% time.seriesID.other]
        
        dat.resurvey1.unassigned = as.data.table(dat.resurvey1)[time_seriesID %in% time.seriesID.unassigned]
        dat.resurvey2.unassigned = as.data.table(dat.resurvey2)[time_seriesID %in% time.seriesID.unassigned]
        
        # calculate species pool richness and linear trends across all
        species.richness.lm = calculate.species.richness.lm(dat.resurvey1, dat.resurvey2)
        results.temp[[as.character(timespan.temp)]][["trajectory_all"]][["EUNIS1_all"]] = 
          data.frame(trajectory = "all",
                     EUNIS1 = "all",
                     decade1 = year.temp, 
                     decade2 = year.temp + timespan.temp,
                     pool_richness1 = length(unique(dat.resurvey1$`sPlot concept unified`)),
                     pool_richness2 = length(unique(dat.resurvey2$`sPlot concept unified`)), 
                     n = length(unique(dat.resurvey1$time_seriesID)),
                     local_richness1 = dat.resurvey1 %>% 
                       group_by(PlotObservationID) %>% 
                       summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                       pull(local_richness) %>% 
                       mean(),
                     local_richness2 = dat.resurvey2 %>% 
                       group_by(PlotObservationID) %>% 
                       summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                       pull(local_richness) %>% 
                       mean(),
                     local_richness_slope = as.numeric(species.richness.lm["slope.year"]))
        
        if(nrow(dat.resurvey1.stable) > 0){
          species.richness.lm = calculate.species.richness.lm(dat.resurvey1.stable, dat.resurvey2.stable)
          results.temp[[as.character(timespan.temp)]][["trajectory_stable"]][["EUNIS1_all"]] = 
            data.frame(trajectory = "stable",
                       EUNIS1 = "all",
                       decade1 = year.temp, 
                       decade2 = year.temp + timespan.temp,
                       pool_richness1 = length(unique(dat.resurvey1.stable$`sPlot concept unified`)),
                       pool_richness2 = length(unique(dat.resurvey2.stable$`sPlot concept unified`)), 
                       n = length(unique(dat.resurvey1.stable$time_seriesID)),
                       local_richness1 = dat.resurvey1.stable %>% 
                         group_by(PlotObservationID) %>% 
                         summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                         pull(local_richness) %>% 
                         mean(),
                       local_richness2 = dat.resurvey2.stable %>% 
                         group_by(PlotObservationID) %>% 
                         summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                         pull(local_richness) %>% 
                         mean(),
                       local_richness_slope = as.numeric(species.richness.lm["slope.year"]))}
        
        if(nrow(dat.resurvey1.succession) > 0){
          species.richness.lm = calculate.species.richness.lm(dat.resurvey1.succession, dat.resurvey2.succession)
          results.temp[[as.character(timespan.temp)]][["trajectory_succession"]][["EUNIS1_all"]] = 
            data.frame(trajectory = "succession",
                       EUNIS1 = "all",
                       decade1 = year.temp, 
                       decade2 = year.temp + timespan.temp,
                       pool_richness1 = length(unique(dat.resurvey1.succession$`sPlot concept unified`)),
                       pool_richness2 = length(unique(dat.resurvey2.succession$`sPlot concept unified`)), 
                       n = length(unique(dat.resurvey1.succession$time_seriesID)),
                       local_richness1 = dat.resurvey1.succession %>% 
                         group_by(PlotObservationID) %>% 
                         summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                         pull(local_richness) %>% 
                         mean(),
                       local_richness2 = dat.resurvey2.succession %>% 
                         group_by(PlotObservationID) %>% 
                         summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                         pull(local_richness) %>% 
                         mean(),
                       local_richness_slope = as.numeric(species.richness.lm["slope.year"]))}
        
        if(nrow(dat.resurvey1.disturbance) > 0){
          species.richness.lm = calculate.species.richness.lm(dat.resurvey1.disturbance, dat.resurvey2.disturbance)
          results.temp[[as.character(timespan.temp)]][["trajectory_disturbance"]][["EUNIS1_all"]] = 
            data.frame(trajectory = "disturbance",
                       EUNIS1 = "all",
                       decade1 = year.temp, 
                       decade2 = year.temp + timespan.temp,
                       pool_richness1 = length(unique(dat.resurvey1.disturbance$`sPlot concept unified`)),
                       pool_richness2 = length(unique(dat.resurvey2.disturbance$`sPlot concept unified`)), 
                       n = length(unique(dat.resurvey1.disturbance$time_seriesID)),
                       local_richness1 = dat.resurvey1.disturbance %>% 
                         group_by(PlotObservationID) %>% 
                         summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                         pull(local_richness) %>% 
                         mean(),
                       local_richness2 = dat.resurvey2.disturbance %>% 
                         group_by(PlotObservationID) %>% 
                         summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                         pull(local_richness) %>% 
                         mean(),
                       local_richness_slope = as.numeric(species.richness.lm["slope.year"]))}
        
        if(nrow(dat.resurvey1.other) > 0){
          species.richness.lm = calculate.species.richness.lm(dat.resurvey1.other, dat.resurvey2.other)
          results.temp[[as.character(timespan.temp)]][["trajectory_other"]][["EUNIS1_all"]] = 
            data.frame(trajectory = "other",
                       EUNIS1 = "all",
                       decade1 = year.temp, 
                       decade2 = year.temp + timespan.temp,
                       pool_richness1 = length(unique(dat.resurvey1.other$`sPlot concept unified`)),
                       pool_richness2 = length(unique(dat.resurvey2.other$`sPlot concept unified`)), 
                       n = length(unique(dat.resurvey1.other$time_seriesID)),
                       local_richness1 = dat.resurvey1.other %>% 
                         group_by(PlotObservationID) %>% 
                         summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                         pull(local_richness) %>% 
                         mean(),
                       local_richness2 = dat.resurvey2.other %>% 
                         group_by(PlotObservationID) %>% 
                         summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                         pull(local_richness) %>% 
                         mean(),
                       local_richness_slope = as.numeric(species.richness.lm["slope.year"]))}
        
        if(nrow(dat.resurvey1.unassigned) > 0){
          species.richness.lm = calculate.species.richness.lm(dat.resurvey1.unassigned, dat.resurvey2.unassigned)
          results.temp[[as.character(timespan.temp)]][["trajectory_unassigned"]][["EUNIS1_all"]] = 
            data.frame(trajectory = "unassigned",
                       EUNIS1 = "all",
                       decade1 = year.temp, 
                       decade2 = year.temp + timespan.temp,
                       pool_richness1 = length(unique(dat.resurvey1.unassigned$`sPlot concept unified`)),
                       pool_richness2 = length(unique(dat.resurvey2.unassigned$`sPlot concept unified`)), 
                       n = length(unique(dat.resurvey1.unassigned$time_seriesID)),
                       local_richness1 = dat.resurvey1.unassigned %>% 
                         group_by(PlotObservationID) %>% 
                         summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                         pull(local_richness) %>% 
                         mean(),
                       local_richness2 = dat.resurvey2.unassigned %>% 
                         group_by(PlotObservationID) %>% 
                         summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                         pull(local_richness) %>% 
                         mean(),
                       local_richness_slope = as.numeric(species.richness.lm["slope.year"]))}
        
        # species pool richness per EUNIS1 class
        if(nrow(dat.resurvey1) > 0){
          for(eunis1.temp in sort(unique(dat.resurvey1$EUNIS1))){
            
            # across all trajectories
            time_seriesIDs.temp = filter(dat.resurvey1, EUNIS1 == eunis1.temp) %>% pull(time_seriesID)
            dat.resurvey1.temp = filter(dat.resurvey1, time_seriesID %in% time_seriesIDs.temp)
            dat.resurvey2.temp = filter(dat.resurvey2, time_seriesID %in% time_seriesIDs.temp)
            
            species.richness.lm = calculate.species.richness.lm(dat.resurvey1.temp, dat.resurvey2.temp)
            results.temp[[as.character(timespan.temp)]][["trajectory_all"]][[eunis1.temp]] = 
              data.frame(trajectory = "all",
                         EUNIS1 = eunis1.temp,
                         decade1 = year.temp, 
                         decade2 = year.temp + timespan.temp,
                         pool_richness1 = length(unique(dat.resurvey1.temp$`sPlot concept unified`)),
                         pool_richness2 = length(unique(dat.resurvey2.temp$`sPlot concept unified`)),
                         n = length(unique(dat.resurvey1.temp$time_seriesID)),
                         local_richness1 = dat.resurvey1.temp %>% 
                           group_by(PlotObservationID) %>% 
                           summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                           pull(local_richness) %>% 
                           mean(),
                         local_richness2 = dat.resurvey2.temp %>% 
                           group_by(PlotObservationID) %>% 
                           summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                           pull(local_richness) %>% 
                           mean(),
                         local_richness_slope = as.numeric(species.richness.lm["slope.year"]))}}
        
        # across "stable" trajectories
        if(nrow(dat.resurvey1.stable) > 0){
          for(eunis1.temp in sort(unique(dat.resurvey1.stable$EUNIS1))){
            time_seriesIDs.temp = filter(dat.resurvey1.stable, EUNIS1 == eunis1.temp) %>% pull(time_seriesID)
            dat.resurvey1.temp = filter(dat.resurvey1.stable, time_seriesID %in% time_seriesIDs.temp)
            dat.resurvey2.temp = filter(dat.resurvey2.stable, time_seriesID %in% time_seriesIDs.temp)
            
            species.richness.lm = calculate.species.richness.lm(dat.resurvey1.temp, dat.resurvey2.temp)
            results.temp[[as.character(timespan.temp)]][["trajectory_stable"]][[eunis1.temp]] = 
              data.frame(trajectory = "stable",
                         EUNIS1 = eunis1.temp,
                         decade1 = year.temp,
                         decade2 = year.temp + timespan.temp,
                         pool_richness1 = length(unique(dat.resurvey1.temp$`sPlot concept unified`)),
                         pool_richness2 = length(unique(dat.resurvey2.temp$`sPlot concept unified`)),
                         n = length(unique(dat.resurvey1.temp$time_seriesID)),
                         local_richness1 = dat.resurvey1.temp %>% 
                           group_by(PlotObservationID) %>% 
                           summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                           pull(local_richness) %>% 
                           mean(),
                         local_richness2 = dat.resurvey2.temp %>% 
                           group_by(PlotObservationID) %>% 
                           summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                           pull(local_richness) %>% 
                           mean(),
                         local_richness_slope = as.numeric(species.richness.lm["slope.year"]))}}
        
        # across "succession" trajectories
        if(nrow(dat.resurvey1.succession) > 0){
          for(eunis1.temp in sort(unique(dat.resurvey1.succession$EUNIS1))){
            time_seriesIDs.temp = filter(dat.resurvey1.succession, EUNIS1 == eunis1.temp) %>% pull(time_seriesID)
            dat.resurvey1.temp = filter(dat.resurvey1.succession, time_seriesID %in% time_seriesIDs.temp)
            dat.resurvey2.temp = filter(dat.resurvey2.succession, time_seriesID %in% time_seriesIDs.temp)
            
            species.richness.lm = calculate.species.richness.lm(dat.resurvey1.temp, dat.resurvey2.temp)
            results.temp[[as.character(timespan.temp)]][["trajectory_succession"]][[eunis1.temp]] = 
              data.frame(trajectory = "succession",
                         EUNIS1 = eunis1.temp,
                         decade1 = year.temp,
                         decade2 = year.temp + timespan.temp,
                         pool_richness1 = length(unique(dat.resurvey1.temp$`sPlot concept unified`)),
                         pool_richness2 = length(unique(dat.resurvey2.temp$`sPlot concept unified`)),
                         n = length(unique(dat.resurvey1.temp$time_seriesID)),
                         local_richness1 = dat.resurvey1.temp %>% 
                           group_by(PlotObservationID) %>% 
                           summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                           pull(local_richness) %>% 
                           mean(),
                         local_richness2 = dat.resurvey2.temp %>% 
                           group_by(PlotObservationID) %>% 
                           summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                           pull(local_richness) %>% 
                           mean(),
                         local_richness_slope = as.numeric(species.richness.lm["slope.year"]))}}
        
        # across "disturbance" trajectories
        if(nrow(dat.resurvey1.disturbance) > 0){
          for(eunis1.temp in sort(unique(dat.resurvey1.disturbance$EUNIS1))){
            time_seriesIDs.temp = filter(dat.resurvey1.disturbance, EUNIS1 == eunis1.temp) %>% pull(time_seriesID)
            dat.resurvey1.temp = filter(dat.resurvey1.disturbance, time_seriesID %in% time_seriesIDs.temp)
            dat.resurvey2.temp = filter(dat.resurvey2.disturbance, time_seriesID %in% time_seriesIDs.temp)
            
            species.richness.lm = calculate.species.richness.lm(dat.resurvey1.temp, dat.resurvey2.temp)
            results.temp[[as.character(timespan.temp)]][["trajectory_disturbance"]][[eunis1.temp]] = 
              data.frame(trajectory = "disturbance",
                         EUNIS1 = eunis1.temp,
                         decade1 = year.temp,
                         decade2 = year.temp + timespan.temp,
                         pool_richness1 = length(unique(dat.resurvey1.temp$`sPlot concept unified`)),
                         pool_richness2 = length(unique(dat.resurvey2.temp$`sPlot concept unified`)),
                         n = length(unique(dat.resurvey1.temp$time_seriesID)),
                         local_richness1 = dat.resurvey1.temp %>% 
                           group_by(PlotObservationID) %>% 
                           summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                           pull(local_richness) %>% 
                           mean(),
                         local_richness2 = dat.resurvey2.temp %>% 
                           group_by(PlotObservationID) %>% 
                           summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                           pull(local_richness) %>% 
                           mean(),
                         local_richness_slope = as.numeric(species.richness.lm["slope.year"]))}}
        
        # across "other" trajectories
        if(nrow(dat.resurvey1.other) > 0){
          for(eunis1.temp in sort(unique(dat.resurvey1.other$EUNIS1))){
            time_seriesIDs.temp = filter(dat.resurvey1.other, EUNIS1 == eunis1.temp) %>% pull(time_seriesID)
            dat.resurvey1.temp = filter(dat.resurvey1.other, time_seriesID %in% time_seriesIDs.temp)
            dat.resurvey2.temp = filter(dat.resurvey2.other, time_seriesID %in% time_seriesIDs.temp)
            
            species.richness.lm = calculate.species.richness.lm(dat.resurvey1.temp, dat.resurvey2.temp)
            results.temp[[as.character(timespan.temp)]][["trajectory_other"]][[eunis1.temp]] = 
              data.frame(trajectory = "other",
                         EUNIS1 = eunis1.temp,
                         decade1 = year.temp,
                         decade2 = year.temp + timespan.temp,
                         pool_richness1 = length(unique(dat.resurvey1.temp$`sPlot concept unified`)),
                         pool_richness2 = length(unique(dat.resurvey2.temp$`sPlot concept unified`)),
                         n = length(unique(dat.resurvey1.temp$time_seriesID)),
                         local_richness1 = dat.resurvey1.temp %>% 
                           group_by(PlotObservationID) %>% 
                           summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                           pull(local_richness) %>% 
                           mean(),
                         local_richness2 = dat.resurvey2.temp %>% 
                           group_by(PlotObservationID) %>% 
                           summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                           pull(local_richness) %>% 
                           mean(),
                         local_richness_slope = as.numeric(species.richness.lm["slope.year"]))}}
        
        # across "unassigned" trajectories
        if(nrow(dat.resurvey1.unassigned) > 0){
          for(eunis1.temp in sort(unique(dat.resurvey1.unassigned$EUNIS1))){
            time_seriesIDs.temp = filter(dat.resurvey1.unassigned, EUNIS1 == eunis1.temp) %>% pull(time_seriesID)
            dat.resurvey1.temp = filter(dat.resurvey1.unassigned, time_seriesID %in% time_seriesIDs.temp)
            dat.resurvey2.temp = filter(dat.resurvey2.unassigned, time_seriesID %in% time_seriesIDs.temp)
            
            species.richness.lm = calculate.species.richness.lm(dat.resurvey1.temp, dat.resurvey2.temp)
            results.temp[[as.character(timespan.temp)]][["trajectory_unassigned"]][[eunis1.temp]] = 
              data.frame(trajectory = "unassigned",
                         EUNIS1 = eunis1.temp,
                         decade1 = year.temp,
                         decade2 = year.temp + timespan.temp,
                         pool_richness1 = length(unique(dat.resurvey1.temp$`sPlot concept unified`)),
                         pool_richness2 = length(unique(dat.resurvey2.temp$`sPlot concept unified`)),
                         n = length(unique(dat.resurvey1.temp$time_seriesID)),
                         local_richness1 = dat.resurvey1.temp %>% 
                           group_by(PlotObservationID) %>% 
                           summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                           pull(local_richness) %>% 
                           mean(),
                         local_richness2 = dat.resurvey2.temp %>% 
                           group_by(PlotObservationID) %>% 
                           summarise(local_richness = length(unique(`sPlot concept unified`))) %>% 
                           pull(local_richness) %>% 
                           mean(),
                         local_richness_slope = as.numeric(species.richness.lm["slope.year"]))}}
        
        print(paste0("done: ", year.temp, " - ", year.temp + timespan.temp))
      }
    }
    return(results.temp)
  }


################################################################################
analyse.and.compile.lms = function(dat.temp, metric.temp){
  
  # must again order by recording_year
  dat.temp = dat.temp[order(dat.temp$year_of_recording),]
  
  # run linear model
  lm.temp = lm(div_value ~ year_of_recording, data = dat.temp)
  
  # compile results
  return(c("metric" = metric.temp,
           "time_seriesID" = dat.temp$time_seriesID[1],
           "intercept" = as.character(lm.temp$coefficients[1]),
           "slope" = as.character(lm.temp$coefficients[2]),
           "initial_div_value" = dat.temp$div_value[1],
           "EUNIS_start" = dat.temp$Expert_system[1],
           "EUNIS_end" = dat.temp$Expert_system[nrow(dat.temp)],
           "year_start" = dat.temp$year_of_recording[1],
           "year_end" = dat.temp$year_of_recording[nrow(dat.temp)],
           "n" = nrow(dat.temp)))}

################################################################################
analyse.and.compile.lms.of.rate.change = function(dat.temp, metric.temp){
  
  # must again sort by year_of_recording
  dat.temp = dat.temp[order(dat.temp$year_of_recording),]
  
  # run linear model
  lm.temp = lm(log(div_value) ~ year_of_recording, data = dat.temp)
  
  # compile results
  return(c("metric" = metric.temp,
           "time_seriesID" = dat.temp$time_seriesID[1],
           "intercept" = as.character(exp(lm.temp$coefficients[1]) - 1),
           "slope" = as.character(exp(lm.temp$coefficients[2]) - 1),
           "initial_div_value" = dat.temp$div_value[1],
           "EUNIS_start" = dat.temp$Expert_system[1],
           "EUNIS_end" = dat.temp$Expert_system[nrow(dat.temp)],
           "year_start" = dat.temp$year_of_recording[1],
           "year_end" = dat.temp$year_of_recording[nrow(dat.temp)],
           "n" = nrow(dat.temp)))}
