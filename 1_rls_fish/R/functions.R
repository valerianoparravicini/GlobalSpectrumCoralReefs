unique_size_sst <- function(rls_fish, rls_predictors) {
  left_join(rls_fish, select(rls_predictors, SurveyID, median_sst_1year)) %>%
    select(Species, size = TL, sst = median_sst_1year) %>%
    drop_na(sst) %>%
    mutate(size = round(size), sst = round(sst)) %>%
    unique() 
}

add_params <- function(sp_size_sst, params) {
  params <- params %>%
    mutate(sst = round(v_m))
  
  linf <- select(params, Species, linf_m) %>% unique()
  
  test <- left_join(sp_size_sst, linf) %>%
    filter(size<linf_m)
  
  inner_join(test, params)
}

run_fishflux <- function(data, cores) {
  parallel::mclapply(1:nrow(data), function(x){
    print(x)
    
    dt <- data[x,] 
    par <- dt %>% select(-species,- size, -Family, - Species, - Genus, - diet, - sst) %>% as.list()
    mod <- fishflux::cnp_model_mcmc(TL = dt$size,
                                    param = par, iter = 1000)
    
    
    extr <- fishflux::extract(mod, par = c("F0c", "F0n", "F0p", "Gc", "Gn", "Gp", "Sc", "Sn", "Sp", 
                                           "Ic", "In", "Ip", "Wc", "Wn", "Wp", "Fc", "Fn", "Fp"))
    extr <- cbind(dt[,1:7], extr) 
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

ind_functions <- function(fluxes_sp_size_sst, troph, sp_size_sst_params) {
  left_join(fluxes_sp_size_sst, troph) %>%
    left_join(unique(select(sp_size_sst_params, species, lwa_m, lwb_m))) %>%
    mutate(biom = lwa_m*(TL^lwb_m),
           herb = p2_m * Ic_median,
           plank = p8_m * Ic_median,
           pisc = p4_m * Ic_median,
           prod = Gc_median,
           exP = Fp_median,
           egP = Wp_median,
           exN = Fn_median,
           egN = Wn_median,
           respC = Fc_median,
           respC_plank = 
             case_when(trophic_guild_predicted == 8 ~ Fc_median,
                                          TRUE ~ 0),
           prod_plank = 
             case_when(trophic_guild_predicted == 8 ~ Gc_median,
                       TRUE ~ 0)
           ) %>%
    mutate(plank = case_when(p8_m < 0.05 ~ 0,
                             TRUE ~ plank),
           herb = case_when(p2_m < 0.05 ~ 0,
                             TRUE ~ herb),
           pisc = case_when(p4_m < 0.05 ~ 0,
                             TRUE ~ pisc)) %>%
    mutate(prod_noplank = prod - prod_plank,
           respC_noplank = respC - respC_plank) %>%
    select(Family, Species, species, TL, sst, biom, herb, plank, pisc, 
           prod, exP, egP, exN, egN, Ic = Ic_median,
           respC, respC_noplank, prod_noplank)  %>%
    mutate(TL = round(TL), sst = round(sst))
}

combine_transects_functions <- function(rls_fish, rls_predictors, rls_meta,  functions_sp_size_sst) {
  data_sst <- rls_predictors %>%
    left_join(select(rls_meta, SiteCode, Lat_Zone)) %>%
    filter(Lat_Zone == "Tropical") %>%
    group_by(SiteCode) %>%
    summarize(sst = round(mean(median_sst_1year, na.rm = T)))
  
  transect <- dplyr::inner_join(rls_fish, data_sst) %>%
    mutate(TL = round(TL)) %>%
    filter(TL>=5) %>%
    as.data.frame() %>%
    dplyr::left_join(as.data.frame(select(functions_sp_size_sst, - Family))) %>%
    mutate(biom2 = biom * Abundance) %>%
    group_by(SurveyID) %>%
    mutate(prop_biom = sum(biom2, na.rm = T)/ sum(Biomass, na.rm = T)) %>%
    filter(prop_biom > 0.7)
  
  length(unique(transect$SurveyID))
  transect
}

summarize_transect <- function(transect) {
  transect %>%
    group_by(SurveyID) %>%
    summarize(
      biom = (sum(biom * Abundance, na.rm = T)) , 
      herb = (sum(herb * Abundance, na.rm = T)) , 
      plank = (sum(plank * Abundance, na.rm = T)) , 
      pisc = (sum(pisc * Abundance, na.rm = T)) , 
      prod = (sum(prod * Abundance, na.rm = T)) , 
      exP = (sum(exP * Abundance, na.rm = T)) , 
      egP = (sum(egP * Abundance, na.rm = T)) , 
      exN = (sum(exN * Abundance, na.rm = T)) , 
      egN = (sum(egN * Abundance, na.rm = T)) , 
      Ic = (sum(Ic * Abundance, na.rm = T)),
      respC = (sum(respC * Abundance, na.rm = T)),
      respC_noplank = (sum(respC_noplank * Abundance, na.rm = T)),
      prod_noplank = (sum(prod_noplank * Abundance, na.rm = T))) %>%
    ungroup() %>%
    mutate(turn = prod/biom,
           exNP = exN/exP,
           outN = exN + egN,
           outP = exP + egP,
           outNP = outN/outP)
}

make_pca <- function(transect_fish_benthic, rls_predictors) {

  # data wrangling
  comb <- transect_fish_benthic %>%
    left_join(select(rls_predictors, SiteCode, Realm, min_sst_5year, max_sst_5year, median_sst_1year, gravtot3)) %>%
    drop_na(gravtot3) %>%
    mutate(sst = median_sst_1year) %>%
    drop_na(sst) %>%
    group_by(Realm) %>%
    mutate(pristine = case_when(gravtot3 < quantile(gravtot3, 0.1) ~ 1, TRUE ~ 0),
           impact = case_when(gravtot3 > quantile(gravtot3, 0.9) ~ 1, TRUE ~ 0)) %>%
    select(- SurveyDate, - SurveyID) %>%
    group_by(SiteCode, Realm, pristine, impact) %>%
    summarize_if(is.numeric, median, na.rm = T) %>%
    ungroup()
    
  
  comb_sel <- comb %>%
    filter(herb>0, plank>0) %>%
    mutate(plank = residuals(lm(log(plank/biom) ~ sst)), 
           herb = residuals(lm(log(herb/biom) ~ sst)), 
           prod = residuals(lm(log(prod/biom) ~ sst)),  
           pisc = residuals(lm(log(pisc/biom) ~ sst)),  
           exP = residuals(lm(log(exP/biom) ~  sst)), 
           exN = residuals(lm(log(exN/biom) ~ sst))) %>%
    mutate(outNP = log(outNP)) 
  
  comb_fish <- comb_sel %>%
    ungroup() %>%
    select(herb, prod, plank, exP, exNP) %>%
    drop_na(herb, prod, plank, exP, exNP)
  
  comb_ben <- comb_sel %>%
    ungroup() %>%
    select(BiomassStand_g, Calc_ghr, BranchSpace_cm3, Rugosity,
           GPP_ghr) %>%
    drop_na()
  
  # pca
  pca_fish <- prcomp(comb_fish, scale = T, center = T) 
  pca_ben <- prcomp(log1p(comb_ben), scale = T, center = T) 
  
  
  pc_fish <- as.data.frame(pca_fish$x) %>%
    cbind(comb_sel)
  rot_fish <- as.data.frame(pca_fish$rotation)
  
  pc_ben <- as.data.frame(pca_ben$x) %>%
    cbind(comb_sel)
  rot_ben <- as.data.frame(pca_ben$rotation)
  
  # plotting
  rot1_ben <- ggplot(rot_ben) +
    geom_vline(aes(xintercept = 0)) +
    geom_hline(aes(yintercept = 0)) +
    geom_point(aes(x = PC1, y = PC2), data = rot_ben) +
    geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), data = rot_ben) +
    geom_label(aes(x = PC1, y = PC2, label = rownames(rot_ben)), size = 3) +
    coord_equal() +
    xlim(c(-1,1)) + ylim(c(-1,1)) +
    theme_classic()

  rot1_fish <- ggplot(rot_fish) +
    geom_vline(aes(xintercept = 0)) +
    geom_hline(aes(yintercept = 0)) +
    geom_point(aes(x = PC1, y = PC2), data = rot_fish) +
    geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), data = rot_fish) +
    geom_label(aes(x = PC1, y = PC2, label = rownames(rot_fish)), size = 3) +
    coord_equal() +
    xlim(c(-1,1)) + ylim(c(-1,1)) +
    theme_classic()

  
  plot <- 
  rot1_fish +     rot1_ben + 

    ggplot(pc_fish, aes(x = PC1, y = PC2)) +
    geom_point(aes(x = PC1, y = PC2),
               data = pc_fish, 
               alpha = 0.4, size = 2) +
    ggforce::geom_mark_hull(aes(x = PC1, y = PC2),
                            data = pc_fish[pc_fish$pristine == 1,],
                            concavity=100,radius=unit(.15,"cm"), color = "blue",
                            expand=unit(.15,"cm"),alpha=.15,size=1)+
    ggforce::geom_mark_hull(aes(x = PC1, y = PC2),
                            data = pc_fish[pc_fish$impact == 1,],
                            concavity=100,radius=unit(.15,"cm"), color = "red",
                            expand=unit(.15,"cm"),alpha=.15,size=1)+
    geom_text(aes(x = 5, y = 4, label = "Pristine"), color = "blue", size = 4) +
    geom_text(aes(x = 5, y = 3, label = "Impacted"), color = "red", size = 4) +
    geom_vline(aes(xintercept = 0)) +
    geom_hline(aes(yintercept = 0))+
    coord_equal() +
    xlim(c(-6,6)) +
    ylim(c(-6,6)) +
    theme_classic() +
    
    
    ggplot(pc_ben, aes(x = PC1, y = PC2)) +
    geom_point(aes(x = PC1, y = PC2),
               data = pc_ben, 
               alpha = 0.4, size = 2) +
    ggforce::geom_mark_hull(aes(x = PC1, y = PC2),
                            data = pc_ben[pc_ben$pristine == 1,],
                            concavity=100,radius=unit(.15,"cm"), color = "blue",
                            expand=unit(.15,"cm"),alpha=.15,size=1)+
    ggforce::geom_mark_hull(aes(x = PC1, y = PC2),
                            data = pc_ben[pc_ben$impact == 1,],
                            concavity=100,radius=unit(.15,"cm"), color = "red",
                            expand=unit(.15,"cm"),alpha=.15,size=1)+
    geom_vline(aes(xintercept = 0)) +
    geom_hline(aes(yintercept = 0))+
    coord_equal() +
    xlim(c(-6,6)) +
    ylim(c(-6,6)) +
    theme_classic() & 
    
    theme(text = element_text(size = 12))
  plot
  ggsave("output/fig1_pca_uglyversion.png", plot, width = 8, height = 8)
  
  ############ per Realm ##############
  plot_realm <- 
    
  ggplot(pc_fish, aes(x = PC1, y = PC2)) +
    geom_point(aes(x = PC1, y = PC2),
               data = pc_fish, 
               alpha = 0.4, size = 2) +
    ggforce::geom_mark_hull(aes(x = PC1, y = PC2),
                            data = pc_fish[pc_fish$pristine == 1,],
                            concavity=100,radius=unit(.15,"cm"), color = "blue",
                            expand=unit(.15,"cm"),alpha=.15,size=1)+
    ggforce::geom_mark_hull(aes(x = PC1, y = PC2),
                            data = pc_fish[pc_fish$impact == 1,],
                            concavity=100,radius=unit(.15,"cm"), color = "red",
                            expand=unit(.15,"cm"),alpha=.15,size=1)+
    geom_vline(aes(xintercept = 0)) +
    geom_hline(aes(yintercept = 0))+
    coord_equal() +
    xlim(c(-6,6)) +
    ylim(c(-6,6)) +
    facet_wrap(~Realm, ncol = 1) +
    labs(title = "Fish") +
    theme_classic() +
    
    ggplot(pc_ben, aes(x = PC1, y = PC2)) +
    geom_point(aes(x = PC1, y = PC2),
               data = pc_ben, 
               alpha = 0.4, size = 2) +
    ggforce::geom_mark_hull(aes(x = PC1, y = PC2),
                            data = pc_ben[pc_ben$pristine == 1,],
                            concavity=100,radius=unit(.15,"cm"), color = "blue",
                            expand=unit(.15,"cm"),alpha=.15,size=1)+
    ggforce::geom_mark_hull(aes(x = PC1, y = PC2),
                            data = pc_ben[pc_ben$impact == 1,],
                            concavity=100,radius=unit(.15,"cm"), color = "red",
                            expand=unit(.15,"cm"),alpha=.15,size=1)+
    geom_vline(aes(xintercept = 0)) +
    geom_hline(aes(yintercept = 0))+
    coord_equal() +
    xlim(c(-6,6)) +
    ylim(c(-6,6)) +
    facet_wrap(~Realm, ncol = 1) +
    labs(title = "Benthos") +
    theme_classic() 
  
  ggsave("output/fig_pca_realm.png", plot_realm, width = 6, height = 12)
  
}

run_models_gravity <- function(transect_fish_benthic, rls_predictors) {
  # data wrangling
  comb <- transect_fish_benthic %>%
    left_join(select(rls_predictors, SiteCode, Realm, min_sst_5year, max_sst_5year, median_sst_1year, gravtot3)) %>%
    drop_na(gravtot3) %>%
    mutate(sst = median_sst_1year) %>%
    drop_na(sst) 
  
  fit_herb <- brm(log(herb) ~ log(gravtot3) + (1|SiteCode), 
                  data = comb[comb$herb>0,], backend = "cmdstanr", cores = 4)
  fit_plank <- brm(log(plank) ~ log(gravtot3) + (1|SiteCode), 
                   data = comb[comb$plank>0,], backend = "cmdstanr", cores = 4)
  fit_exP <- brm(log(exP) ~ log(gravtot3) + (1|SiteCode), 
                 data = comb, backend = "cmdstanr", cores = 4)
  fit_exNP <- brm(log(exNP) ~ log(gravtot3) + (1|SiteCode), 
                  data = comb, backend = "cmdstanr", cores = 4)
  fit_prod <- brm(log(prod) ~ log(gravtot3) + (1|SiteCode), 
                  data = comb, backend = "cmdstanr", cores = 4)
  
  fit_calc <- brm(log(Calc_ghr) ~ log(gravtot3) + (1|SiteCode), 
                  data = comb[comb$Calc_ghr > 0,], backend = "cmdstanr", cores = 4)
  fit_space <- brm(log(BranchSpace_cm3) ~ log(gravtot3) + (1|SiteCode), 
                   data = comb[comb$BranchSpace_cm3 > 0,], backend = "cmdstanr", cores = 4)
  fit_sbio <- brm(log(BiomassStand_g) ~ log(gravtot3) + (1|SiteCode), 
                  data = comb, backend = "cmdstanr", cores = 4)
  fit_rug <- brm(log(Rugosity) ~ log(gravtot3) + (1|SiteCode), 
                 data = comb, backend = "cmdstanr", cores = 4)
  fit_gpp <- brm(log(GPP_ghr) ~ log(gravtot3) + (1|SiteCode), 
                 data = comb, backend = "cmdstanr", cores = 4)
  
  return(list(fit_prod, fit_exP, fit_exNP, fit_plank, fit_herb,
              fit_calc, fit_space, fit_sbio, fit_rug, fit_gpp))
  # summary(fit_herb)
  # summary(fit_plank)
  # summary(fit_exP)
  # summary(fit_exNP)
  # summary(fit_prod)
  # 
  # summary(fit_calc)
  # summary(fit_space)
  # summary(fit_sbio)
  # summary(fit_rug)
  # summary(fit_gpp)
  
}

run_models_coralcover <- function(transect_fish_benthic, rls_predictors) {
  # data wrangling
  comb <- transect_fish_benthic %>%
    left_join(select(rls_predictors, SiteCode, Realm, min_sst_5year, max_sst_5year, median_sst_1year, gravtot3)) %>%
    drop_na(gravtot3) %>%
    mutate(sst = median_sst_1year) %>%
    drop_na(sst) 
  
  fit_herb <- brm(log(herb) ~ sqrt(coral) + (1|SiteCode), 
                  data = comb[comb$herb>0,], backend = "cmdstanr", cores = 4)
  fit_plank <- brm(log(plank) ~ sqrt(coral) + (1|SiteCode), 
                   data = comb[comb$plank>0,], backend = "cmdstanr", cores = 4)
  fit_exP <- brm(log(exP) ~ sqrt(coral) + (1|SiteCode), 
                 data = comb, backend = "cmdstanr", cores = 4)
  fit_exNP <- brm(log(exNP) ~ sqrt(coral) + (1|SiteCode), 
                  data = comb, backend = "cmdstanr", cores = 4)
  fit_prod <- brm(log(prod) ~ sqrt(coral) + (1|SiteCode), 
                  data = comb, backend = "cmdstanr", cores = 4)
  
  fit_calc <- brm(log(Calc_ghr) ~ sqrt(coral) + (1|SiteCode), 
                  data = comb[comb$Calc_ghr > 0,], backend = "cmdstanr", cores = 4)
  fit_space <- brm(log1p(BranchSpace_cm3) ~ sqrt(coral) + (1|SiteCode), 
                   data = comb, backend = "cmdstanr", cores = 4)
  fit_sbio <- brm(log(BiomassStand_g) ~ sqrt(coral) + (1|SiteCode), 
                  data = comb, backend = "cmdstanr", cores = 4)
  fit_rug <- brm(log(Rugosity) ~ sqrt(coral) + (1|SiteCode), 
                 data = comb, backend = "cmdstanr", cores = 4)
  fit_gpp <- brm(log(GPP_ghr) ~ sqrt(coral) + (1|SiteCode), 
                 data = comb, backend = "cmdstanr", cores = 4)
  
  return(list(fit_prod, fit_exP, fit_exNP, fit_plank, fit_herb,
              fit_calc, fit_space, fit_sbio, fit_rug, fit_gpp))

  
}

extract_slopes <- function(mod_list) {
  lapply(mod_list, function(x){fixef(x)[2,]}) %>% plyr::ldply() %>%
    mutate(fun = c("Production", "P excretion", "N:P flux", "Planktivory", "Herbivory",
                   "Calcification", "Branch space", "Standing biomass", "Rugosity", "GPP"))
    
}


make_plot2 <- function(transect_fish_benthic, rls_predictors, models_grav, models_cover, slopes_grav, slopes_cover) {
  
  comb <- transect_fish_benthic %>%
    left_join(select(rls_predictors, SiteCode, Realm, min_sst_5year, max_sst_5year, median_sst_1year, gravtot3)) %>%
    drop_na(gravtot3) %>%
    mutate(sst = median_sst_1year) %>%
    drop_na(sst) 
  
  pred1 <- lapply(models_grav, function(x){
    fitted(x, re_formula = NA)
  })
  pred2 <- lapply(models_cover, function(x){
    fitted(x, re_formula = NA)
  })
  
  dens <- lapply(models_grav, function(x){
    pr <- fitted(x, re_formula = NA, 
                 newdata = data.frame(gravtot3 = c(quantile(comb$gravtot3, 0.05),
                                                   quantile(comb$gravtot3, 0.95))), 
                 summary = F)
    ggplot() +
      geom_halfeyeh(aes( x = exp(pr[,1])), 
                    color = "blue", fill = "blue", alpha = 0.5) +
      geom_halfeyeh(aes( x = exp(pr[,2])), 
                    color = "red", fill = "red", alpha = 0.5) +
      scale_y_continuous(breaks = c(0,1)) +
      theme_classic() +
      theme(
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'), 
        title = element_blank()
      )
  })
  
  ggplot(comb) +
    geom_point(aes(x = log(gravtot3), y = log(prod)), 
                   alpha = 0.5, size = 0.5) +
    geom_ribbon(aes(x = log(gravtot3), ymin = pred1[[1]][,3], ymax = pred1[[1]][,4]), 
                alpha = 0.5)  +
    geom_line(aes(x = log(gravtot3), y = pred1[[1]][,1]))  +
    labs(x = "Human gravity", y = "log(Production)", 
         title = paste0("slope = ", round(slopes_grav[1, 1], 3), " (",  
                        round(slopes_grav[1, 3], 3), " ; ", 
                        round(slopes_grav[1, 4], 3), ")")) +
    theme_classic() +
    theme(legend.position = "top") +
    
    inset_element(dens[[1]], 0.7, 0.7, 1,1) +
    
    ggplot() +
    geom_point(aes(x = log(gravtot3), y = log(Calc_ghr)),
               alpha = 0.5, size = 0.5, data = comb) +
    geom_ribbon(aes(x = models_grav[[6]]$data$`log(gravtot3)`, ymin = pred1[[6]][,3], ymax = pred1[[6]][,4]), 
                alpha = 0.5)  +
    geom_line(aes(x = models_grav[[6]]$data$`log(gravtot3)`, y = pred1[[6]][,1]))  +
    labs(x = "Human gravity", y = "log(Calcification)",
         title = paste0("slope = ", round(slopes_grav[6, 1], 3), " (",  
                        round(slopes_grav[6, 3], 3), " ; ", 
                        round(slopes_grav[6, 4], 3), ")")) +
    theme_classic() +
    theme(legend.position = "none") +
    
    inset_element(dens[[6]], 0.7, 0.7, 1,1) +
    
    ggplot(comb) +
    geom_point(aes(x = log(gravtot3), y = log(exP)),
               alpha = 0.5, size = 0.5) +
    geom_ribbon(aes(x = log(gravtot3), ymin = pred1[[2]][,3], ymax = pred1[[2]][,4]), 
                alpha = 0.5)  +
    geom_line(aes(x = log(gravtot3), y = pred1[[2]][,1]))  +
    labs(x = "Human gravity", y = "log(P excretion)",
         title = paste0("slope = ", round(slopes_grav[2, 1], 3), " (",  
                        round(slopes_grav[2, 3], 3), " ; ", 
                        round(slopes_grav[2, 4], 3), ")")) +
    theme_classic() +
    theme(legend.position = "none") +
    
    inset_element(dens[[2]], 0.7, 0.7, 1,1) +
    
    ggplot(comb) +
    geom_point(aes(x = log(gravtot3), y = log(BranchSpace_cm3)),
               alpha = 0.5, size = 0.5) +
    geom_ribbon(aes(x = log(gravtot3), ymin = pred1[[7]][,3], ymax = pred1[[7]][,4]), 
                alpha = 0.5)  +
    geom_line(aes(x = log(gravtot3), y = pred1[[7]][,1]))  +
    labs(x = "Human gravity", y = "log(Space creation)",
         title = paste0("slope = ", round(slopes_grav[7, 1], 3), " (",  
                        round(slopes_grav[7, 3], 3), " ; ", 
                        round(slopes_grav[7, 4], 3), ")")) +
    theme_classic() +
    theme(legend.position = "none") +
    
    inset_element(dens[[7]], 0.7, 0.7, 1,1) +
    
    ggplot(comb) +
    geom_point(aes(x = log(gravtot3), y = log(exNP)),
               alpha = 0.5, size = 0.5) +
    geom_ribbon(aes(x = log(gravtot3), ymin = pred1[[3]][,3], ymax = pred1[[3]][,4]), 
                alpha = 0.5)  +
    geom_line(aes(x = log(gravtot3), y = pred1[[3]][,1]))  +
    labs(x = "Human gravity", y = "log(N:P flux)",
         title = paste0("slope = ", round(slopes_grav[3, 1], 3), " (",  
                        round(slopes_grav[3, 3], 3), " ; ", 
                        round(slopes_grav[3, 4], 3), ")")) +
    theme_classic() +
    theme(legend.position = "none") +
    
    inset_element(dens[[3]], 0.7, 0.7, 1,1) +
    
    ggplot(comb) +
    geom_point(aes(x = log(gravtot3), y = log(BiomassStand_g)),
               alpha = 0.5, size = 0.5) +
    geom_ribbon(aes(x = log(gravtot3), ymin = pred1[[8]][,3], ymax = pred1[[8]][,4]), 
                alpha = 0.5)  +
    geom_line(aes(x = log(gravtot3), y = pred1[[8]][,1]))  +
    labs(x = "Human gravity", y = "log(Standing biomass)",
         title = paste0("slope = ", round(slopes_grav[8, 1], 3), " (",  
                        round(slopes_grav[8, 3], 3), " ; ", 
                        round(slopes_grav[8, 4], 3), ")")) +
    theme_classic() +
    theme(legend.position = "none") +
    
    inset_element(dens[[8]], 0.7, 0.7, 1,1) +
    
    ggplot() +
    geom_point(aes(x = log(gravtot3), y = log(plank)),
               alpha = 0.5, size = 0.5, data = comb) +
    geom_ribbon(aes(x =models_grav[[4]]$data$`log(gravtot3)`, ymin = pred1[[4]][,3], ymax = pred1[[4]][,4]), 
                alpha = 0.5)  +
    geom_line(aes(x = models_grav[[4]]$data$`log(gravtot3)`, y = pred1[[4]][,1]))  +
    labs(x = "Human gravity", y = "log(Planktivory)",
         title = paste0("slope = ", round(slopes_grav[4, 1], 3), " (",  
                        round(slopes_grav[4, 3], 3), " ; ", 
                        round(slopes_grav[4, 4], 3), ")")) +
    theme_classic() +
    theme(legend.position = "none") +
    
    inset_element(dens[[4]], 0.7, 0.7, 1,1) +
    
    ggplot(comb) +
    geom_point(aes(x = log(gravtot3), y = log(Rugosity)),
               alpha = 0.5, size = 0.5) +
    geom_ribbon(aes(x = log(gravtot3), ymin = pred1[[9]][,3], ymax = pred1[[9]][,4]), 
                alpha = 0.5)  +
    geom_line(aes(x = log(gravtot3), y = pred1[[9]][,1]))  +
    labs(x = "Human gravity", y = "log(Rugosity)",
         title = paste0("slope = ", round(slopes_grav[9, 1], 3), " (",  
                        round(slopes_grav[9, 3], 3), " ; ", 
                        round(slopes_grav[9, 4], 3), ")")) +
    theme_classic() +
    theme(legend.position = "none") +
    
    inset_element(dens[[9]], 0.7, 0.7, 1,1) +
    
   
    ggplot() +
    geom_point(aes(x = log(gravtot3), y = log(herb)),
               alpha = 0.5, size = 0.5, data = comb) +
    geom_ribbon(aes(x = models_grav[[5]]$data$`log(gravtot3)`, ymin = pred1[[5]][,3], ymax = pred1[[5]][,4]), 
                alpha = 0.5)  +
    geom_line(aes(x = models_grav[[5]]$data$`log(gravtot3)`, y = pred1[[5]][,1]))  +
    labs(x = "Human gravity", y = "log(Herbivory)",
         title = paste0("slope = ", round(slopes_grav[5, 1], 3), " (",  
                        round(slopes_grav[5, 3], 3), " ; ", 
                        round(slopes_grav[5, 4], 3), ")")) +
    theme_classic() +
    theme(legend.position = "none") +
    
    inset_element(dens[[5]], 0.7, 0.7, 1,1) +
    
    ggplot(comb) +
    geom_point(aes(x = log(gravtot3), y = log(GPP_ghr)),
               alpha = 0.5, size = 0.5) +
    geom_ribbon(aes(x = log(gravtot3), ymin = pred1[[10]][,3], ymax = pred1[[10]][,4]), 
                alpha = 0.5)  +
    geom_line(aes(x = log(gravtot3), y = pred1[[10]][,1]))  +
    labs(x = "Human gravity", y = "log(GPP)",
         title = paste0("slope = ", round(slopes_grav[10, 1], 3), " (",  
                        round(slopes_grav[10, 3], 3), " ; ", 
                        round(slopes_grav[10, 4], 3), ")")) +
    theme_classic() +
    theme(legend.position = "none") +
    
    inset_element(dens[[10]], 0.7, 0.7, 1,1) +
    
    plot_layout(ncol = 2)
    
  ggsave("output/fig2_gravity_uglyversion.png", width = 12, height = 18)
  
  
  ggplot(comb) +
    geom_point(aes(x = sqrt(coral), y = log(prod)), 
               alpha = 0.5, size = 0.5) +
    geom_ribbon(aes(x = sqrt(coral), ymin = pred2[[1]][,3], ymax = pred2[[1]][,4]), 
                alpha = 0.5)  +
    geom_line(aes(x = sqrt(coral), y = pred2[[1]][,1]))  +
    labs(x = "sqrt(Coral cover)", y = "log(Production)", 
         title = paste0("slope = ", round(slopes_cover[1, 1], 3), " (",  
                        round(slopes_cover[1, 3], 3), " ; ", 
                        round(slopes_cover[1, 4], 3), ")")) +
    theme_classic() +
    theme(legend.position = "top") +
    
    ggplot() +
    geom_point(aes(x = sqrt(coral), y = log(Calc_ghr)),
               alpha = 0.5, size = 0.5, data = comb) +
    geom_ribbon(aes(x = models_cover[[6]]$data$`sqrt(coral)`, ymin = pred2[[6]][,3], ymax = pred2[[6]][,4]), 
                alpha = 0.5)  +
    geom_line(aes(x = models_cover[[6]]$data$`sqrt(coral)`, y = pred2[[6]][,1]))  +
    labs(x = "sqrt(Coral cover)", y = "log(Calcification)",
         title = paste0("slope = ", round(slopes_cover[6, 1], 3), " (",  
                        round(slopes_cover[6, 3], 3), " ; ", 
                        round(slopes_cover[6, 4], 3), ")")) +
    theme_classic() +
    theme(legend.position = "none") +
    
    ggplot(comb) +
    geom_point(aes(x = sqrt(coral), y = log(exP)),
               alpha = 0.5, size = 0.5) +
    geom_ribbon(aes(x = sqrt(coral), ymin = pred2[[2]][,3], ymax = pred2[[2]][,4]), 
                alpha = 0.5)  +
    geom_line(aes(x = sqrt(coral), y = pred2[[2]][,1]))  +
    labs(x = "sqrt(Coral cover)", y = "log(P excretion)",
         title = paste0("slope = ", round(slopes_cover[2, 1], 3), " (",  
                        round(slopes_cover[2, 3], 3), " ; ", 
                        round(slopes_cover[2, 4], 3), ")")) +
    theme_classic() +
    theme(legend.position = "none") +
    
    ggplot(comb) +
    geom_point(aes(x = sqrt(coral), y = log(BranchSpace_cm3)),
               alpha = 0.5, size = 0.5) +
    geom_ribbon(aes(x = sqrt(coral), ymin = pred2[[7]][,3], ymax = pred2[[7]][,4]), 
                alpha = 0.5)  +
    geom_line(aes(x = sqrt(coral), y = pred2[[7]][,1]))  +
    labs(x = "sqrt(Coral cover)", y = "log(Space creation)",
         title = paste0("slope = ", round(slopes_cover[7, 1], 3), " (",  
                        round(slopes_cover[7, 3], 3), " ; ", 
                        round(slopes_cover[7, 4], 3), ")")) +
    theme_classic() +
    theme(legend.position = "none") +
    
    ggplot(comb) +
    geom_point(aes(x = sqrt(coral), y = log(exNP)),
               alpha = 0.5, size = 0.5) +
    geom_ribbon(aes(x = sqrt(coral), ymin = pred2[[3]][,3], ymax = pred2[[3]][,4]), 
                alpha = 0.5)  +
    geom_line(aes(x = sqrt(coral), y = pred2[[3]][,1]))  +
    labs(x = "sqrt(Coral cover)", y = "log(N:P flux)",
         title = paste0("slope = ", round(slopes_cover[3, 1], 3), " (",  
                        round(slopes_cover[3, 3], 3), " ; ", 
                        round(slopes_cover[3, 4], 3), ")")) +
    theme_classic() +
    theme(legend.position = "none") +
    
    ggplot(comb) +
    geom_point(aes(x = sqrt(coral), y = log(BiomassStand_g)),
               alpha = 0.5, size = 0.5) +
    geom_ribbon(aes(x = sqrt(coral), ymin = pred2[[8]][,3], ymax = pred2[[8]][,4]), 
                alpha = 0.5)  +
    geom_line(aes(x = sqrt(coral), y = pred2[[8]][,1]))  +
    labs(x = "sqrt(Coral cover)", y = "log(Standing biomass)",
         title = paste0("slope = ", round(slopes_cover[8, 1], 3), " (",  
                        round(slopes_cover[8, 3], 3), " ; ", 
                        round(slopes_cover[8, 4], 3), ")")) +
    theme_classic() +
    theme(legend.position = "none") +
    
    ggplot() +
    geom_point(aes(x = sqrt(coral), y = log(plank)),
               alpha = 0.5, size = 0.5, data = comb) +
    geom_ribbon(aes(x =models_cover[[4]]$data$`sqrt(coral)`, ymin = pred2[[4]][,3], ymax = pred2[[4]][,4]), 
                alpha = 0.5)  +
    geom_line(aes(x = models_cover[[4]]$data$`sqrt(coral)`, y = pred2[[4]][,1]))  +
    labs(x = "sqrt(Coral cover)", y = "log(Planktivory)",
         title = paste0("slope = ", round(slopes_cover[4, 1], 3), " (",  
                        round(slopes_cover[4, 3], 3), " ; ", 
                        round(slopes_cover[4, 4], 3), ")")) +
    theme_classic() +
    theme(legend.position = "none") +
    
    ggplot(comb) +
    geom_point(aes(x = sqrt(coral), y = log(Rugosity)),
               alpha = 0.5, size = 0.5) +
    geom_ribbon(aes(x = sqrt(coral), ymin = pred2[[9]][,3], ymax = pred2[[9]][,4]), 
                alpha = 0.5)  +
    geom_line(aes(x = sqrt(coral), y = pred2[[9]][,1]))  +
    labs(x = "sqrt(Coral cover)", y = "log(Rugosity)",
         title = paste0("slope = ", round(slopes_cover[9, 1], 3), " (",  
                        round(slopes_cover[9, 3], 3), " ; ", 
                        round(slopes_cover[9, 4], 3), ")")) +
    theme_classic() +
    theme(legend.position = "none") +
    
    
    ggplot() +
    geom_point(aes(x = sqrt(coral), y = log(herb)),
               alpha = 0.5, size = 0.5, data = comb) +
    geom_ribbon(aes(x = models_cover[[5]]$data$`sqrt(coral)`, ymin = pred2[[5]][,3], ymax = pred2[[5]][,4]), 
                alpha = 0.5)  +
    geom_line(aes(x = models_cover[[5]]$data$`sqrt(coral)`, y = pred2[[5]][,1]))  +
    labs(x = "sqrt(Coral cover)", y = "log(Herbivory)",
         title = paste0("slope = ", round(slopes_cover[5, 1], 3), " (",  
                        round(slopes_cover[5, 3], 3), " ; ", 
                        round(slopes_cover[5, 4], 3), ")")) +
    theme_classic() +
    theme(legend.position = "none") +
    
    ggplot(comb) +
    geom_point(aes(x = sqrt(coral), y = log(GPP_ghr)),
               alpha = 0.5, size = 0.5) +
    geom_ribbon(aes(x = sqrt(coral), ymin = pred2[[10]][,3], ymax = pred2[[10]][,4]), 
                alpha = 0.5)  +
    geom_line(aes(x = sqrt(coral), y = pred2[[10]][,1]))  +
    labs(x = "sqrt(Coral cover)", y = "log(GPP)",
         title = paste0("slope = ", round(slopes_cover[10, 1], 3), " (",  
                        round(slopes_cover[10, 3], 3), " ; ", 
                        round(slopes_cover[10, 4], 3), ")")) +
    theme_classic() +
    theme(legend.position = "none") +
    plot_layout(ncol = 2)
  
  ggsave("output/fig2_coral_uglyversion.png", width = 8, height = 12)
  
}

plot_percent_change <- function() {
  comb <- transect_fish_benthic %>%
    left_join(select(rls_predictors, SiteCode, Realm, min_sst_5year, max_sst_5year, median_sst_1year, gravtot3)) %>%
    drop_na(gravtot3) %>%
    mutate(sst = median_sst_1year) %>%
    drop_na(sst) 
  
  plts <- lapply(models_grav, function(x){
    pr <- fitted(x, re_formula = NA, 
                 newdata = data.frame(gravtot3 = c(quantile(comb$gravtot3, 0.05),
                                                   quantile(comb$gravtot3, 0.95))), 
                 summary = F)
    ggplot() +
      geom_halfeyeh(aes( x = exp(pr[,1])), 
                    color = "blue", fill = "blue", alpha = 0.5) +
      geom_halfeyeh(aes( x = exp(pr[,2])), 
                    color = "red", fill = "red", alpha = 0.5) +
      theme_classic() +
      labs(x = "", y = "")
  })
  
  
  
}
