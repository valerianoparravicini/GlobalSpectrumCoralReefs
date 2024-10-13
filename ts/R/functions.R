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
  

  
  ####### models #######

  fit_prod <- brm(log(prod) ~ time  + (time|location)+ (time|site)+ (time|site_zone),
                 data = data[data$prod > 0, ], backend = "cmdstanr", cores = 4, iter = 4000)
  
  fit_exP <- brm(log(exP) ~ time  + (time|location) + (time|site)+ (time|site_zone),
                 data = data, backend = "cmdstanr", cores = 4, iter = 4000)
  
  fit_exNP <- brm(log(exNP) ~ time  + (time|location)+ (time|site)+ (time|site_zone),
                  data = data, backend = "cmdstanr", cores = 4, iter = 4000)
  
  fit_plank <- brm(log(plank) ~ time  + (time|location)+ (time|site)+ (time|site_zone),
                   data = data[data$plank>0,], backend = "cmdstanr", cores = 4, iter = 4000)
  
  fit_herb <- brm(log(herb) ~ time  + (time|location)+ (time|site)+ (time|site_zone),
                   data = data[data$herb>0,], backend = "cmdstanr", cores = 4, iter = 4000)
  
  list(fit_prod, fit_exP, fit_exNP, fit_plank, fit_herb)
  
  # test <-  cbind(fit_plank$data,
  #   fitted(fit_plank, re_formula = "log(plank) ~ time  + (time|location)"))
  # 
  # nd <-  data.frame(time = c(0,1))
  # test2 <- cbind(nd, fitted(fit_plank, re_formula = "log(plank) ~ time", newdata = nd))
  # 
  # ggplot(test) +
  #   geom_line(aes(x = time, y = (Estimate), group = location ),
  #             alpha = 0.5) +
  #   geom_ribbon(aes(x = time, ymin = `Q2.5`, ymax = `Q97.5`),
  #               data = test2, size = 2, fill = "blue", alpha = 0.2) +
  #   geom_line(aes(x = time, y = Estimate), data = test2, size = 2, color = "blue") +
  #   theme_bw()
  #   
  # 

  
}

plot_ts01 <- function(ts_models_01) {
  
  f <- c(" Production", "P excretion", "N:P Flux", "Planktivory", "Herbivory")
  
  plot_list <- lapply(1:5, function(i){
    
    x <- ts_models_01[[i]]
    name <- f[i]

    nd0 <-  select(x$data, location, site, site_zone) %>% unique() %>%
      mutate(time = 0) 
  
    nd1 <-  select(x$data, location, site, site_zone) %>% unique() %>%
      mutate(time = 1) 
    
    pred0 <- fitted(x, newdata = nd0, summary = F)
    pred1 <- fitted(x, newdata = nd1, summary = F)
    pred <-exp(pred1)/exp(pred0)
    
    nd1$pred <- apply(pred, 2, median)
    nd1$lb <- apply(pred, 2, quantile, 0.025)
    nd1$ub <- apply(pred, 2, quantile, 0.975)
    
    cat <- nd1 %>%
      mutate(cat = case_when(lb>1 ~ "increase",
                             ub<1 ~ "decrease",
                             TRUE ~ "rest")) %>%
      select(location, site, site_zone, cat)
    
    nd0$pred <- 1
    nd0$lb <- 1
    nd0$ub <- 1
    
    ndpred <- rbind(nd0, nd1) %>%
      left_join(cat)
    
    b <- fixef(x, pars = "time", summary = F)
    
    pdens <- ggplot() +
      geom_halfeyeh(aes(x = b), alpha = 0.5) +
      theme_classic() +
      scale_y_continuous(breaks = c(0,1)) +
      geom_vline(aes(xintercept = 0), linetype = 2) +
      labs(x = "Global slope", y = "Density") +
      theme(
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')
      ) 
    
    
    plot <- ggplot() +
      geom_line(aes(x = time, y = (pred), group = paste0(location, site, site_zone)),
                alpha = 0.25, data = ndpred[ndpred$cat == "rest",]) +
      geom_line(aes(x = time, y = (pred), group = paste0(location, site,site_zone)),
                alpha = 0.5, data = ndpred[ndpred$cat == "increase",],
                color = "blue") +
      geom_line(aes(x = time, y = (pred), group = paste0(location, site, site_zone)),
                alpha = 0.5, data = ndpred[ndpred$cat == "decrease",],
                color = "red") +
      scale_y_continuous(trans = "log10") +
      theme_classic() +
      labs(x = "Time", y = name) +
      inset_element(pdens, 0,0.7, 0.3, 1)
      
      return(plot)
    
  })
  
  plot_list[[1]]  +
    plot_list[[2]]  +
    plot_list[[3]]  +
    plot_list[[4]]  +
    plot_list[[5]]  +
    plot_layout(ncol = 1) 

ggsave("output/time01_overview.png", width = 12, height = 18)


}

