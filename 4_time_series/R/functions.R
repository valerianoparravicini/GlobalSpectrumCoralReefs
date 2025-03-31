clean_fish <- function(fish_raw) {
  
  fish <- fish_raw %>%
    janitor::clean_names() %>%
    filter(!is.na(size)) %>%
    mutate(size = round(size)) %>%
    dplyr::group_by(dataset_id, location, site, zone) %>%
    dplyr::mutate(nyear = length(unique(year))) %>%
    filter(nyear > 2) %>%
    mutate(site = str_replace_all(site, "[^[:alnum:]]", "")) %>%
    mutate(species = gsub(" ", "_", species)) %>%
    mutate(location = gsub(" ", "", location),
           site = gsub(" ", "", site)) %>%
    ungroup() 
  
  # add a dataset
  mpa <- read_csv("data/moorea_uvc_mpa.csv") %>%
    separate(taxon, into = c("genus", "sp"), remove = FALSE) %>%
    mutate(dataset_id = "Moorea_MPA",
           replicate = rep_id,
           area = "Pacific",
           country = "France",
           archipelago = "Society Islands",
           location = "Moorea",
           site = site_name,
           zone = reef_zone,
           latitude = lat,
           longitude = long,
           date = lubridate::as_date(paste(year, month, day, sep = "-")),
           method = "uvc_100",
           observer = "criobe",
           species = gsub(" ", "_", taxon),
           density = abundance
    ) %>%
    dplyr::group_by(dataset_id, location, site, zone) %>%
    dplyr::mutate(nyear = length(unique(year))) 
  
  moo_mpa <- mpa[, colnames(mpa) %in% colnames(fish)]
  
  fish %>%
    full_join(moo_mpa) %>%
    filter(size > 4) 
  
}



clean_bent <- function(benthic_raw) {
  benthic_raw %>% janitor::clean_names() %>%
    filter(category == "coral") %>%
    mutate(location = gsub(" ", "", location),
           site = gsub(" ", "", site)) %>%
    group_by(dataset_id, country, location, site, latitude, longitude, replicate, quadrat, zone, year, depth) %>%
    summarise(coral_cover = sum(cover)) %>%
    mutate(site = str_replace_all(site, "[^[:alnum:]]", "")) %>%
    mutate(location = case_when(location == "MontebelloBarrowIslands" ~ "Montebellos",
                                TRUE ~ location)) %>%
    ungroup()
}



  

unique_spsize <- function(fish) {
  fish %>% 
    filter(size > 0) %>%
    mutate(size = round(size)) %>%
    select(species, size) %>%
    ungroup() %>%
    unique() 
}


sst_time_series <- function(data) {
  data <- janitor::clean_names(data)
  sst_dataset <- "ncdcOisst21Agg_LonPM180"
  sst_info <- rerddap::info(
    sst_dataset, url = "https://coastwatch.pfeg.noaa.gov/erddap/"
  )
  # Produce a dataset with daily values for each site.
  # Info from `info(sst_dataset)`:
  #   var is "sst", dims of interest are latitude, longitude and time
  all_coords <- data %>%
    dplyr::distinct(year, latitude, longitude) %>%
    split(f = ~ latitude + longitude, drop = TRUE) %>%
    purrr::map_dfr(function(x) {
      data.frame(latitude = x$latitude[1], longitude = x$longitude[1],
                 year = min(x$year):max(x$year))
    })
  all_ssts <- vector(mode = "list", length = nrow(all_coords))
  global_expansion <- c(0.05, -0.05)
  for (j in seq_along(all_ssts)) {
    latitude <- global_expansion + all_coords$latitude[j]
    longitude <- global_expansion + all_coords$longitude[j]
    yr_j <- all_coords$year[j]
    day_1 <- paste0(yr_j, "-01-01")
    day_2 <- paste0(yr_j, ifelse(yr_j == 2021, "-11-11", "-12-31"))
    dates <- c(day_1, day_2)
    sst_j <- try(rerddap::griddap(
      sst_info, latitude = latitude, longitude = longitude, time = dates,
      fields = "sst"
    ), silent = TRUE)
    n <- 1
    mod_expansion <- global_expansion + c(0.05, -0.05)
    while (inherits(sst_j, "try-error")) {
      if (n > 10) {
        break
      }
      latitude <- mod_expansion + all_coords$latitude[j]
      longitude <- mod_expansion + all_coords$longitude[j]
      message("row ", j, "; expanding box by 0.05 deg; attempt ", n)
      sst_j <- try(rerddap::griddap(
        sst_info, latitude = latitude, longitude = longitude, time = dates,
        fields = "sst"
      ), silent = TRUE)
      n <- n + 1
      mod_expansion <- mod_expansion + c(0.05, -0.05)
    }
    if (n > 10 && inherits(sst_j, "try-error")) {
      next
    }
    all_ssts[[j]] <- sst_j %>%
      `[[`("data") %>%
      dplyr::group_by(time) %>%
      dplyr::summarise(sst = mean(sst, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::summarise(year = all_coords$year[j],  sst = mean(sst))
    attr(all_ssts[[j]], "info") <- all_coords[j, ]
  }
  # checks
  # which(sapply(all_ssts, is.null))
  # sum(sapply(all_ssts, is.null))
  names(all_ssts) <- sapply(all_ssts, function(x) {
    attr(x, "info") %>%
      unlist %>%
      paste0(collapse = "_:_")
  }, USE.NAMES = FALSE)
  all_ssts
}

extract_yr <- function(x) {
  format(extract_time(x), "%Y")
}


unique_spsize_sst <- function(fish, sst_data) {
  fish %>% 
    filter(size > 0) %>%
    mutate(size = round(size)) %>%
    left_join(sst_data) %>%
    mutate(sst = round(sst)) %>%
    select(species, size, v_m = sst) %>%
    ungroup() %>%
    unique() %>%
    drop_na()
}


add_params <- function(sp_size, params) {
  
  params <- dplyr::select(params, -size)
  linf <- select(params, species, linf_m) %>% unique()
  
  test <- left_join(sp_size, linf) %>%
    mutate(big = size>linf_m*2) %>%
    mutate(size = case_when(big ~ linf_m, TRUE ~ size)) %>%
    select(-big)
  
  left_join(test, params) %>%
    ungroup() %>%
    drop_na(Family) # only keep combinations for which we have parameters (~99%)
}

run_fishflux <- function(data, cores) {
  parallel::mclapply(1:nrow(data), function(x){
    print(x)
    
    dt <- data[x,] 
    par <- dt %>% select(-c(species, size, Family, Genus, Species, diet)) %>% as.list()
    mod <- fishflux::cnp_model_mcmc(TL = dt$size,
                                    param = par, iter = 1000)
    
    
    extr <- fishflux::extract(mod, par = c("F0c", "F0n", "F0p", "Gc", "Gn", "Gp", "Sc", "Sn", "Sp", 
                                           "Ic", "In", "Ip", "Wc", "Wn", "Wp", "Fc", "Fn", "Fp"))
    extr <- cbind(select(dt, c(species, size, Family, Genus, Species, diet)), extr) 
    lim <- fishflux::limitation(mod, plot = FALSE)
    extr$limitation <-first(lim[lim$prop_lim == max(lim$prop_lim), "nutrient"])
    
    return(extr)
  }, mc.cores = cores) %>% plyr::ldply()
}

run_fishflux_segmented <- function(data, cores) {
  
  start <- seq(1, nrow(data), 500)
  end <- c(start[-1] - 1, nrow(data)) 
  rows <- data.frame(start, end)
  
  lapply(1:nrow(rows), function(i){
    print(paste("#####################", rows[i,1], ":", rows[i,2], "######################"))
    run_fishflux(data[rows[i,1]:rows[i,2],], cores = cores)
  } ) %>% plyr::ldply()
  
}

ind_functions <- function(fluxes_sp_size, troph, sp_size_params) {
  left_join(fluxes_sp_size, troph) %>%
    left_join(unique(select(sp_size_params, species, lwa_m, lwb_m))) %>%
    mutate(biom = lwa_m*(TL^lwb_m),
           herb = p2_m * Ic_median,
           plank = p8_m * Ic_median,
           pisc = p4_m * Ic_median,
           prod = Gc_median,
           exP = Fp_median,
           egP = Wp_median,
           exN = Fn_median,
           egN = Wn_median) %>%
    mutate(plank = case_when(p8_m < 0.1 ~ 0,
                             TRUE ~ plank),
           herb = case_when(p2_m < 0.1 ~ 0,
                            TRUE ~ herb),
           pisc = case_when(p4_m < 0.1 ~ 0,
                            TRUE ~ pisc)) %>%
    select(Family, Species, species, TL, biom, herb, plank, pisc, 
           prod, exP, egP, exN, egN, Ic = Ic_median)  %>%
    mutate(TL = round(TL))
}

combine_transects_functions <- function(fish, functions_sp_size) {
  
  functions <- functions_sp_size %>%
    group_by(species, TL) %>%
    summarise_if(is.numeric, median) %>%
    mutate(size = TL)
  
  transect <- dplyr::inner_join(fish, functions) %>%
    filter(size >= 5) %>%
    as.data.frame() 
  
  transect
}

summarize_transect <- function(transect, coral) {
  
  cor <- select(coral, location, site, zone, year, coral_cover) %>%
    group_by(location, year, site, zone) %>%
    summarize(coral_cover = mean(coral_cover))
  
  transect %>%
    group_by(dataset_id, area, country, archipelago, location, site, replicate, zone, 
             latitude, longitude, depth, year, nyear) %>%
    summarize(
      biom = sum(biom * density, na.rm = T) , 
      herb = sum(herb * density, na.rm = T) , 
      plank = sum(plank * density, na.rm = T) , 
      pisc = sum(pisc * density, na.rm = T) , 
      prod = sum(prod * density, na.rm = T) , 
      exP = sum(exP * density, na.rm = T) , 
      egP = sum(egP * density, na.rm = T) , 
      exN = sum(exN * density, na.rm = T) , 
      egN = sum(egN * density, na.rm = T) , 
      Ic = sum(Ic * density, na.rm = T)) %>%
    ungroup() %>%
    mutate(turn = prod/biom,
           egNP = egN/egP,
           exNP = exN/exP,
           outN = exN + egN,
           outP = exP + egP,
           outNP = outN/outP) %>%
    left_join(cor)
}

fit_ts_models_fish <- function(data) {
  
  #data <- summary_fish_functions 
  
  data[data$location == "Kon\xe9" , "location"] <- "Kone"
  
  data <- data %>%
    group_by(location, dataset_id, site, zone) %>%
    drop_na(location) %>%
    mutate(time = (year - min(year))/(max(year) - min(year)))  %>%
    mutate(loc_site = paste(location, site, sep = "_")) %>%
    mutate(site_zone = paste(site, zone, sep = "_"))
  
}
  

  
# benthos

clean_benthos_ts <- function(ts_benthos_case_study) {

# rock and rubble act more as turf?
ts_benthos_case_study$taxon_original <- ts_benthos_case_study$taxon
ts_benthos_case_study$taxon <- ifelse(ts_benthos_case_study$taxon=="Rock" | ts_benthos_case_study$taxon=="coral rubble", "turf algae", as.character(ts_benthos_case_study$taxon))

# macroalgae acts similar to fleshy? 
ts_benthos_case_study$taxon <- ifelse(ts_benthos_case_study$taxon=="canopy forming macroalgae" | ts_benthos_case_study$taxon=="understory macroalgae", "fleshy algae", as.character(ts_benthos_case_study$taxon))

ts_benthos_case_study$N <- ceiling(ts_benthos_case_study$cover) # % cover = number of "individuals"
sites <- unique(paste(ts_benthos_case_study$Year, ts_benthos_case_study$Site, ts_benthos_case_study$Replicate, sep = "_"))
ts_benthos_case_study$SurveyID = paste(ts_benthos_case_study$Year, ts_benthos_case_study$Site, ts_benthos_case_study$Replicate, sep = "_")
ntot   <-length(unique(sites))
coords <- expand.grid(c(1:10), c(1:10))  # for plotting

# Matching RLS previous dataset
ts_benthos_case_study$taxon[ts_benthos_case_study$taxon == "Hard_Coral_Branching"] = "Branching"
ts_benthos_case_study$taxon[ts_benthos_case_study$taxon == "Hard_Coral_Tabular"] = "Tabular"
ts_benthos_case_study$taxon[ts_benthos_case_study$taxon == "Hard_Coral_Encrusting"] = "Encrusting"
ts_benthos_case_study$taxon[ts_benthos_case_study$taxon %in% c("Hard_Coral_Massive", "Hard_Coral_Submassive", "Massive")] = "Hemispherical"
ts_benthos_case_study$taxon[ts_benthos_case_study$taxon %in% c("Alcyonacea", "soft")] = "ahermatypic coral"

ts_benthos_case_study

}#eo clean_benthos_ts 

  
check_XY_difference <- function(x, y) { c(FALSE, abs(diff(x)) == 1 | abs(diff(y)) == 1) }


sim_colony_size_ts <- function(ts_benthos_case_study) {
  
  sites <- unique(paste(ts_benthos_case_study$Year, ts_benthos_case_study$Site, ts_benthos_case_study$Replicate, sep = "_"))
  ntot   <-length(unique(sites))
  coords <- expand.grid(c(1:10), c(1:10)) 
  
  indivs <- parallel::mclapply(sites, function(i) {
    
    set.seed(10)
    site <- ts_benthos_case_study[ts_benthos_case_study$SurveyID==i  ,]
    inds <- rep(site$taxon, site$N)
    diam_cm <- 10
    xy <- coords[sample(nrow(coords)),][1:length(inds),] # randomise for plotting
    data.frame(SurveyID=i, taxon = inds, diam_cm, x=xy$Var1, y=xy$Var2)
    
  }, mc.cores=50)
  
  indivs <- do.call(rbind, indivs)
  

  # Cumulated size
  
  surveys <- length(unique(indivs$SurveyID))
  surveys_dataset <- indivs %>% group_by(SurveyID) %>% group_split()
  
  surveys_size_arranged <- vector(mode = "list", length = surveys)
  
  for (i in 1:surveys) {
    surveys_size_arranged[[i]] <- surveys_dataset[[i]] %>%
      group_by(taxon) %>%
      mutate(group_id = cumsum(!check_XY_difference(x, y))) %>%
      group_by(taxon, group_id) %>%
      summarise(SurveyID = first(SurveyID),
                X = mean(x),
                Y = mean(y),
                diam_cm = sum(diam_cm)) 
    }
  
  indivs = surveys_size_arranged %>% bind_rows()
  
  indivs
  
}#eo sim_colony_size_ts


hemisphere<-function(d, tb, ...){
  r <- (d/2)
  SA<-(2*pi*r^2)
  TB<-SA*tb
  VOL<-((2/3)*pi*r^3)
  FULL <- VOL # No branch space
  data.frame(SA, TB, VOL, FULL)}

flat<-function(d, th, tb, ...){
  r <- (d/2)
  SA<-pi*(r^2)
  TB<-SA*tb
  VOL<-SA*th
  FULL <- VOL # No branch space
  data.frame(SA, TB, VOL, FULL)}

branches<-function(d, bh, br, bpa, tb, ...){ 
  r <- (d/2)
  SA<-(pi*(r^2)*(bpa*(2*pi*br*bh)+(pi*(br^2))))
  TB<-SA*tb
  VOL<-(pi*(r^2)*(bpa*(pi*(br^2)*bh)))
  FULL <- (pi*((d/2)^2)) * bh
  data.frame(SA,TB, VOL, FULL)}

cone<-function(d, bh, th, tb, ...){
  r <- d/2
  #SA <- pi*r*(r+sqrt((bh^2)+(r^2)))
  SA <- pi * r * sqrt((r^2)+(bh^2)) #pi*r*l
  TB<-SA*tb
  VOL<-th*(SA) 
  FULL <- (pi*(r^2)) * bh
  data.frame(SA,TB, VOL, FULL)}



geometry_change_1_year_ts <- function(traits, indivs, geo, rls_info) {

growth <- subset(traits, trait_name=="Growth rate")
growth <- aggregate(value~rls_group, growth, mean) 


indivs$growth <- growth$value[match(indivs$taxon, growth$rls_group)]
indivs$diam_cm_next <- indivs$diam_cm + (indivs$growth /10 )
head(indivs, 30)

indivs$taxon <- ifelse(indivs$taxon=="Rock and Rubble" | indivs$taxon=="Rubble", "Rock and Rubble", as.character(indivs$taxon))
indivs$taxon <- ifelse(indivs$taxon=="Macroalgae", "fleshy algae", as.character(indivs$taxon))


# set formula
indivs$fun <- rls_info$geometry[match(indivs$taxon, rls_info$taxon)]
use <- indivs[complete.cases(indivs), ] # removes all without trait data 
head(use, 30)

use <- use %>% data.frame()

t0 <- NULL
t1 <- NULL
ntot <- nrow(use)

for (i in 1:ntot){
  funct <- use[i, "fun"]
  gp <- use[i,"taxon"]
  pars <- geo[geo$rls==gp, c("br","bh", "bpa", "th", "tb")]
  pars$d <- use[i,"diam_cm"]
  t0 <- rbind(t0, unlist(do.call(Map, c(f= as.name(funct), pars))))
  pars$d <- use[i,"diam_cm_next"]
  t1 <- rbind(t1, unlist(do.call(Map, c(f= as.name(funct), pars))))
  print(round(which(c(1:ntot)==i)/ntot*100, 2))
}
colnames(t1) <- paste(colnames(t1), "_next", sep="")
morphs <- cbind(use, t0, t1)
head(morphs)

# Rugosity corrected from Husband et al. 2021 (https://doi.org/10.1007/s00338-022-02247-6)
morphs$Rugosity[morphs$taxon %in% c("Halimeda", "fleshy algae", "turf algae")] <- 1
morphs$Rugosity[morphs$taxon %in% c("Encrusting", "coralline algae")] <- 1.6
morphs$Rugosity[morphs$taxon == "Hemispherical"] <- 1.8
morphs$Rugosity[morphs$taxon == "Tabular"] <- 2.0
morphs$Rugosity[morphs$taxon == "Corymbose"] <- 2.3
morphs$Rugosity[morphs$taxon == "Branching"] <- 3.0
morphs$Rugosity[morphs$taxon == "Laminar"] <- 3.5

# Rugosity per SurveyID
morphs_Rugo <- morphs %>% mutate(Rugosity_reef = (Rugosity * diam_cm)) %>% group_by(SurveyID) %>% 
  summarise(Rugosity_reef = mean(Rugosity_reef)/10)

hist(morphs_Rugo$Rugosity_reef)

# Merge Rugosity to previous data
morphs <- merge(morphs, morphs_Rugo) %>% dplyr::select(-Rugosity) %>% rename(Rugosity = Rugosity_reef)

# merge additional traits with morphs data 
trait_avs <- aggregate(value~rls_group+trait_name, traits, mean) # or sample?

return(list(morphs,  trait_avs, morphs_Rugo))

}#eo geometry_change_1_year_ts




algae_temp_factor <- function(toptprime, erprime, ei, temp_c) {
  topt <- 288.15 + plogis(toptprime) * 30
  er <- ei * plogis(erprime)
  boltz <- 8.62e-5
  ei * exp(
    er * ((temp_c + 273.15) - topt) / ((temp_c + 273.15) * boltz * topt)
  ) /
    (ei - er * (1 - exp(
      ei * ((temp_c + 273.15) - topt) / ((temp_c + 273.15) * boltz * topt))
    ))
}

coral_temp_factor <- function(g_max = 1, rho, temp_shift_c, ln_sigma, temp_c) {
  g_max * exp( -1 * exp(rho * (temp_c - temp_shift_c)) -
                 exp(ln_sigma) * (temp_c - temp_shift_c)^2
  )
}



compute_functions <- function(morphs_traits, predictors_ts, rls_info, ts_benthos_case_study = ts_benthos_case_study_clean) {
  
  morphs <- morphs_traits[[1]]
  trait_avs <- morphs_traits[[2]]
  morphs_Rugo <- morphs_traits[[3]] 
  
  ######################## (1) carbonate production & storage
  
  # change in volume and skeletal density
  sd_dat <- subset(trait_avs, trait_name=="Skeletal density")
  morphs$skel_den <- sd_dat$value[match(morphs$taxon, sd_dat$rls_group)]
  morphs$Storage_kg <- (morphs$VOL*morphs$skel_den )/1000 # (cm3 * gcm3 = g )
  morphs$Accretion_kgyr <- ((morphs$VOL_next-morphs$VOL)*morphs$skel_den)/1000 
  morphs$Storage_kg[is.na(morphs$Storage_kg)] <- 0  # others contribute 0
  morphs$Accretion_kgyr[is.na(morphs$Accretion_kgyr)] <- 0 # others contribute 0
  
  ######################## (2) branch space & rugosity
  
  morphs$BranchSpace_cm3<- (morphs$FULL - morphs$VOL) 
  remove <- c("fleshy algae", "turf algae", "Halimeda") # don't contribute?
  morphs$BranchSpace_cm3[morphs$taxon %in% remove ]<-NA
  morphs$BranchSpace_cm3[is.na(morphs$BranchSpace_cm3)] <- 0 # others contribute 0
  
  ######################## (3) photosynthesis & calcification
  
  gpp_dat <- subset(trait_avs, trait_name=="Gross photosynthesis")
  morphs$gpp <- gpp_dat$value[match(morphs$taxon, gpp_dat$rls_group)] #ugcm2hr
  morphs$GPP_ghr <- ( morphs$SA * morphs$gpp) / (1000^2)
  
  calc_dat <- subset(trait_avs, trait_name=="Calcification rate")
  morphs$calc <- calc_dat$value[match(morphs$taxon, calc_dat$rls_group)] #ugcm2hr
  morphs$Calc_ghr <- ( morphs$SA * morphs$calc) / (1000^2)
  morphs$Calc_ghr[is.na(morphs$Calc_ghr)] <- 0 # others contribute 0
  head(morphs)
  
  ######################## (4) organic production & storage
  
  # corals by SA algae by vol
  morphs$calc_status <- rls_info$calcification[match(morphs$taxon, rls_info$taxon)]
  morphs$OrgGrowth_gyr <- (morphs$TB_next - morphs$TB) / 1000 # corals
  morphs$OrgGrowth_gyr <- ifelse(morphs$calc_status=="non_calcifying", 
                                 (((morphs$VOL_next - morphs$VOL)*0.07) - 15) / 1000, morphs$OrgGrowth_gyr)
  # based on https://doi.org/10.1515/BOT.2002.063
  
  morphs$BiomassStand_g <- (morphs$TB) / 1000
  morphs$BiomassStand_g  <- ifelse(morphs$calc_status=="non_calcifying", 
                                   ((morphs$VOL*0.07) - 15) / 1000, morphs$BiomassStand_g) 
  
  # based on https://doi.org/10.1515/BOT.2002.063
  
  pred <- predictors_ts %>% 
    dplyr::select(Year, Location, mean_sst) %>% 
    mutate(merge = paste(Location, Year, sep = "_"))
  
  dataset_sites_merging <- data.frame(Site = ts_benthos_case_study$Site, Location = ts_benthos_case_study$dataset) %>% distinct()
  morphs$Site = gsub("^\\d+_(.*)_\\d+$", "\\1", morphs$SurveyID)
  morphs$Year <- as.numeric(substr(morphs$SurveyID, 1, 4))
  
  morphs <- morphs %>% left_join(dataset_sites_merging) #%>% left_join(pred)
  morphs$merge <- paste(morphs$Location, morphs$Year, sep = "_")
  morphs <- morphs %>% left_join(pred, by = "merge") %>% dplyr::select(-merge)
  
  morphs$coral_sst_factor <- coral_temp_factor( g_max = 1, rho = 1.022, temp_shift_c = 30.197, ln_sigma = -4.518,    temp_c = morphs$mean_sst)
  
  morphs$algae_sst_factor <- algae_temp_factor(
    toptprime = 0.21, erprime = -1.31, ei = 3, temp_c = morphs$mean_sst  )
  
  algae <- c("coralline algae", "fleshy algae", "Halimeda", "turf algae")
  morphs$coral_sst_factor[morphs$taxon %in% algae] <- 1
  morphs$algae_sst_factor[!morphs$taxon %in% algae] <- 1
  rates <- c("OrgGrowth_gyr", "GPP_ghr", "Accretion_kgyr", "Calc_ghr")
  
  morphs[,rates] <- morphs[,rates] * morphs$coral_sst_factor * morphs$algae_sst_factor
  
  functions <- c("Storage_kg", "Accretion_kgyr", "GPP_ghr", "Calc_ghr", "BranchSpace_cm3", "OrgGrowth_gyr", "BiomassStand_g")
  
  fun_taxa <- aggregate(.~SurveyID + taxon, morphs[,c("SurveyID", "taxon", functions)], sum)
  funs <- merge(fun_taxa, morphs_Rugo, by = "SurveyID") %>% rename(Rugosity = Rugosity_reef)
  
  label_decomposition <- str_split(funs$SurveyID, fixed("_"))
  
  label_decomposition <- data.frame(do.call(rbind, label_decomposition))
  label_decomposition[,1] <- as.numeric(label_decomposition[,1]) 
  label_decomposition[,3] <- as.numeric(label_decomposition[,3]) 
  
  
  funs$Location <- label_decomposition[,2] 
  funs$Year <- label_decomposition[,1] 
  
  funs = funs %>% rename(Site = Location) %>% left_join(dataset_sites_merging)
  
  benthos_final = funs %>% dplyr::select(-c(SurveyID, taxon)) %>% group_by(Location, Site, Year) %>% 
    summarise(Storage_avg = mean(Storage_kg), Storage_sd = sd(Storage_kg), 
              Accretion_kgyr_avg = mean(Accretion_kgyr), Accretion_kgyr_sd = sd(Accretion_kgyr), 
              GPP_ghr_avg = mean(GPP_ghr), GPP_ghr_sd = sd(GPP_ghr),
              Calc_ghr_avg = mean(Calc_ghr), Calc_ghr_sd = sd(Calc_ghr),
              BranchSpace_cm3_avg = mean(BranchSpace_cm3), BranchSpace_cm3_sd = sd(BranchSpace_cm3),
              OrgGrowth_gyr_avg = mean(OrgGrowth_gyr), OrgGrowth_gyr_sd = sd(OrgGrowth_gyr),
              BiomassStand_g_avg = mean(BiomassStand_g), BiomassStand_g_sd = sd(BiomassStand_g),
              Rugosity_avg = mean(Rugosity), Rugosity_sd = sd(Rugosity)) 
  
  benthos_final
  
}


create_ts_functions_dataset <- function(fish, benthos) {

#fish <- read_csv("output/ts_transect_fish_functions.csv")
fish$location[fish$location %in% c("Cousin", "Mahe", "Praslin", "Ste")] = "Seychelles"

fish = fish %>% dplyr::select(c(location, site, year, herb, plank, pisc, prod, exP, exN, turn)) %>% group_by(location, site, year) %>% 
  summarise(herb_avg = mean(herb), plank_avg = mean(plank), pisc_avg = mean(pisc), prod_avg = mean(prod), exP_avg = mean(exP), exN_avg = mean(exN), turn_avg = mean(turn))

benthos = benthos %>% rename(site = Site, location = Location, year = Year)
benthos$site <- gsub(" " , "", benthos$site )
fish$site <- gsub(" " , "", fish$site)

ts_functions_dataset <- na.omit(merge(benthos, fish, by = c("year", "location", "site"), all = T))

ts_functions_dataset 

}


#load("../3_rls_analysis/output/pca_functions.RData")

predict_pca <- function(ts_functions_dataset) {
  
  ts_functions_dataset <-  ts_functions_dataset %>% rename(herb = herb_avg, plank = plank_avg, prod = prod_avg, pisc = pisc_avg, exN = exN_avg, exP = exP_avg, turn = turn_avg) 
  ts_functions_dataset <- ts_functions_dataset %>% rename(GPP_ghr = GPP_ghr_avg, OrgGrowth_gyr = OrgGrowth_gyr_avg, Calc_ghr = Calc_ghr_avg, Storage_kg = Storage_avg, Accretion_kgyr = Accretion_kgyr_avg, Rugosity = Rugosity_avg, BranchSpace_cm3 = BranchSpace_cm3_avg) 
  
  fun <- c("herb", "plank", "prod","pisc",  "exN", "exP", "turn", "GPP_ghr", "OrgGrowth_gyr",  "Calc_ghr", "Storage_kg", "Accretion_kgyr",  "Rugosity", "BranchSpace_cm3")
  all <- ts_functions_dataset[, fun]
  
  log.all <- log(all+1)
  
  load("../3_rls_analysis/output/pca_functions.RData")
  
  ts_pca <- data.frame(predict(pca_functions, newdata = log.all))
  
  ts_pca$year <- ts_functions_dataset$year 
  ts_pca$location <- ts_functions_dataset$location
  ts_pca$site <- ts_functions_dataset$site
  
  ts_pca
  
}


