

prepare_PC_data <- function(dat, pred) {


################################# BENTHOS

IDs <- dat$SurveyID

env <- c("SurveyID", "mean_chl_5year", "mean_sss_5year", "mean_npp_5year", "mean_sst_5year", "mean_pH_1year_5year")

environment <- lapply(IDs, function(x) {
  
  na.omit(pred[pred$SurveyID == x, which(colnames(pred) %in% env)])[1,]  
})

env_df <- do.call(rbind, environment) 
env_df <- na.omit(env_df)

dat <- dat[dat$SurveyID %in% env_df$SurveyID,]
pred <- pred[pred$SurveyID %in% env_df$SurveyID,]

pred$No.take.multizoned[pred$No.take.multizoned == "No take"] <- "No take"
pred$No.take.multizoned[pred$No.take.multizoned == "No take "] <- "No take"
pred$No.take.multizoned[pred$No.take.multizoned == "No Take"] <- "No take"
pred$No.take.multizoned[pred$No.take.multizoned == "No take multizoned"] <- "No take"
pred$No.take.multizoned[pred$No.take.multizoned == "Restricted take multizoned"] <- "Restricted take"
pred$No.take.multizoned[pred$No.take.multizoned == "Restricted take"] <- "Restricted take"
pred$No.take.multizoned[pred$No.take.multizoned == "Fishing"] <- "Fishing"
pred$No.take.multizoned[is.na(pred$No.take.multizoned)] <- "Not protected"

IDs <- dat$SurveyID

DHW <- lapply(IDs, function(x) {
  
  c(na.omit(pred[pred$SurveyID == x,]$mean_DHW_5year)[1], na.omit(pred[pred$SurveyID == x,]$mean_DHW_1year[1]))  
})

names(DHW) <- IDs


DHW_df <- do.call(rbind, DHW) 


sites <- sapply(IDs, function(x) {
  
  pred[pred$SurveyID == x,]$SiteCode[1]
  
})

gravtot <- sapply(IDs, function(x) {
  
  pred[pred$SurveyID == x,]$gravtot2[1]
  
})

mpa <- sapply(IDs, function(x) {
  
  pred[pred$SurveyID == x,]$No.take.multizoned[1]
  
})

pc <- prcomp(env_df[,2:6],
             center = TRUE,
             scale. = TRUE)

dat$DHW5 <- DHW_df[,1]
dat$DHW1 <- DHW_df[,2]
dat$gravtot <- gravtot
dat$SiteCode <- sites
dat$mpa <- mpa
dat$env1 <- pc$x[,1]
dat$env2 <- pc$x[,2]

dat

}

###########  PC models env 

run_models_env <- function(dat_mod) {

library(brms)
  
pr <- get_prior(PC1 ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = dat_mod)
pr$prior[which(pr$class == "b")] <- "normal(0,1)"

fit_pc1_env <- brm(PC1 ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = dat_mod, backend = "cmdstanr", cores = 4, threads = 10, iter = 5000, seed = 10, prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))

pr <- get_prior(PC2 ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = dat_mod)
pr$prior[which(pr$class == "b")] <- "normal(0,1)"

fit_pc2_env <- brm(PC2 ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = dat_mod, backend = "cmdstanr", cores = 4, threads = 10, iter = 5000, seed = 10, prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))

pr <- get_prior(PC3 ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = dat_mod)
pr$prior[which(pr$class == "b")] <- "normal(0,1)"

fit_pc3_env <- brm(PC3 ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = dat_mod, backend = "cmdstanr", cores = 4, threads = 10, iter = 5000, seed = 10, prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))

pr <- get_prior(PC4 ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = dat_mod)
pr$prior[which(pr$class == "b")] <- "normal(0,1)"

fit_pc4_env <- brm(PC4 ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = dat_mod, backend = "cmdstanr", cores = 4, threads = 10, iter = 5000, seed = 10, prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))

fit_pc_env <- list(fit_pc1_env,fit_pc2_env, fit_pc3_env, fit_pc4_env)


}


###########  PC models no_env 


run_models_noenv <- function(dat_mod) {
    
    library(brms)
  
    pr <- get_prior(PC1 ~ DHW1  + gravtot + mpa + (1|SiteCode), data = dat_mod)
    pr$prior[which(pr$class == "b")] <- "normal(0,1)"
    
    fit_pc1_noenv <- brm(PC1 ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = dat_mod, backend = "cmdstanr", cores = 4, threads = 10, iter = 5000, seed = 10, prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
   
    pr <- get_prior(PC2 ~ DHW1  + gravtot + mpa + (1|SiteCode), data = dat_mod)
    pr$prior[which(pr$class == "b")] <- "normal(0,1)"
    
    fit_pc2_noenv <- brm(PC2 ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = dat_mod, backend = "cmdstanr", cores = 4, threads = 10, iter = 5000, seed = 10, prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
    
    pr <- get_prior(PC3 ~ DHW1  + gravtot + mpa + (1|SiteCode), data = dat_mod)
    pr$prior[which(pr$class == "b")] <- "normal(0,1)"
    
    fit_pc3_noenv <- brm(PC3 ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = dat_mod, backend = "cmdstanr", cores = 4, threads = 10, iter = 5000, seed = 10, prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
    
    pr <- get_prior(PC4 ~ DHW1  + gravtot + mpa + (1|SiteCode), data = dat_mod)
    pr$prior[which(pr$class == "b")] <- "normal(0,1)"
    
    fit_pc4_noenv <- brm(PC4 ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = dat_mod, backend = "cmdstanr", cores = 4, threads = 10, iter = 5000, seed = 10, prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
    
fit_pc_noenv <- list(fit_pc1_noenv,fit_pc2_noenv, fit_pc3_noenv, fit_pc4_noenv)

fit_pc_noenv
  
}


######## prepare benthic models

prepare_benthic_models <- function(pred, benthos) {
 
################################# BENTHOS

IDsB <- benthos$SurveyID

env <- c("SurveyID", "mean_chl_5year", "mean_sss_5year", "mean_npp_5year", "mean_sst_5year", "mean_pH_1year_5year")
environment <- lapply(IDsB, function(x) {
  
  na.omit(pred[pred$SurveyID == x, which(colnames(pred) %in% env)])[1,]  
})

env_df <- do.call(rbind, environment) 
env_df <- na.omit(env_df)

benthos <- benthos[benthos$SurveyID %in% env_df$SurveyID,]
#fish <- fish[fish$SurveyID %in% env_df$SurveyID,]
pred <- pred[pred$SurveyID %in% env_df$SurveyID,]

pred$No.take.multizoned[pred$No.take.multizoned == "No take"] <- "No take"
pred$No.take.multizoned[pred$No.take.multizoned == "No take "] <- "No take"
pred$No.take.multizoned[pred$No.take.multizoned == "No Take"] <- "No take"
pred$No.take.multizoned[pred$No.take.multizoned == "No take multizoned"] <- "No take"
pred$No.take.multizoned[pred$No.take.multizoned == "Restricted take multizoned"] <- "Restricted take"
pred$No.take.multizoned[pred$No.take.multizoned == "Restricted take"] <- "Restricted take"
pred$No.take.multizoned[pred$No.take.multizoned == "Fishing"] <- "Fishing"
pred$No.take.multizoned[is.na(pred$No.take.multizoned)] <- "Not protected"

IDsB <- benthos$SurveyID


DHWB <- lapply(IDsB, function(x) {
  
  c(na.omit(pred[pred$SurveyID == x,]$mean_DHW_5year)[1], na.omit(pred[pred$SurveyID == x,]$mean_DHW_1year[1]))  
  
  
})

names(DHWB) <- IDsB


DHWB_df <- do.call(rbind, DHWB) 


sites <- sapply(IDsB, function(x) {
  
  pred[pred$SurveyID == x,]$SiteCode[1]
  
})

gravtot <- sapply(IDsB, function(x) {
  
  pred[pred$SurveyID == x,]$gravtot2[1]
  
})

mpa <- sapply(IDsB, function(x) {
  
  pred[pred$SurveyID == x,]$No.take.multizoned[1]
  
})

pc <- prcomp(env_df[,2:6],
             center = TRUE,
             scale. = TRUE)

benthos$DHW5 <- DHWB_df[,1]
benthos$DHW1 <- DHWB_df[,2]
benthos$gravtot <- gravtot
benthos$SiteCode <- sites
benthos$mpa <- mpa
benthos$env1 <- pc$x[,1]
benthos$env2 <- pc$x[,2]

benthos

}





##################################### FISH


prepare_fish_models <- function(pred, fish) {
IDsF <- fish$SurveyID

env <- c("SurveyID", "mean_chl_5year", "mean_sss_5year", "mean_npp_5year", "mean_sst_5year", "mean_pH_1year_5year")
environment <- lapply(IDsF, function(x) {
  
  na.omit(pred[pred$SurveyID == x, which(colnames(pred) %in% env)])[1,]  
})

env_df <- do.call(rbind, environment) 
env_df <- na.omit(env_df)

fish <- fish[fish$SurveyID %in% env_df$SurveyID,]
pred <- pred[pred$SurveyID %in% env_df$SurveyID,]



pred$No.take.multizoned[pred$No.take.multizoned == "No take"] <- "No take"
pred$No.take.multizoned[pred$No.take.multizoned == "No take "] <- "No take"
pred$No.take.multizoned[pred$No.take.multizoned == "No Take"] <- "No take"
pred$No.take.multizoned[pred$No.take.multizoned == "No take multizoned"] <- "No take"
pred$No.take.multizoned[pred$No.take.multizoned == "Restricted take multizoned"] <- "Restricted take"
pred$No.take.multizoned[pred$No.take.multizoned == "Restricted take"] <- "Restricted take"
pred$No.take.multizoned[pred$No.take.multizoned == "Fishing"] <- "Fishing"
pred$No.take.multizoned[is.na(pred$No.take.multizoned)] <- "Not protected"

IDsF <- fish$SurveyID

DHWF <- lapply(IDsF, function(x) {
  
  c(na.omit(pred[pred$SurveyID == x,]$mean_DHW_5year)[1], na.omit(pred[pred$SurveyID == x,]$mean_DHW_1year[1]))  
  
  
})

names(DHWF) <- IDsF

DHWF_df <- do.call(rbind, DHWF) 


sites <- sapply(IDsF, function(x) {
  
  pred[pred$SurveyID == x,]$SiteCode[1]
  
})

gravtot <- sapply(IDsF, function(x) {
  
  pred[pred$SurveyID == x,]$gravtot2[1]
  
})

mpa <- sapply(IDsF, function(x) {
  
  pred[pred$SurveyID == x,]$No.take.multizoned[1]
  
})

environment <- lapply(IDsF, function(x) {
  
  na.omit(pred[pred$SurveyID == x, which(colnames(pred) %in% env)])[1,]  
})

env_df <- do.call(rbind, environment) 
env_df <- na.omit(env_df)

pc <- prcomp(env_df[,2:6],
             center = TRUE,
             scale. = TRUE)

fish$DHW5 <- DHWF_df[,1]
fish$DHW1 <- DHWF_df[,2]
fish$SiteCode <- sites
fish$gravtot <- gravtot
fish$mpa <- mpa
fish$env1 <- pc$x[,1]
fish$env2 <- pc$x[,2]

fish
}


ind_benthic_models_env <- function(benthos) {

###########  IND benthic models env
library(brms)

pr <- get_prior(log1p(Calc_ghr) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = benthos)
pr$prior[which(pr$class == "b")] <- "normal(0,1)"

#fit env

#control = list(adapt_delta = 0.9, max_treedepth = 12),

fit_cc <- brm(log1p(coral) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 10, iter = 5000, seed = 10,  prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
fit_cal <- brm(log1p(Calc_ghr) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 10, iter = 5000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
fit_sto <- brm(log1p(Storage_kg) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 10, iter = 5000, seed = 10,  prior = p, control = list(adapt_delta = 0.9, max_treedepth = 12))
fit_rug <- brm(log1p(Rugosity) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 10, iter = 5000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
fit_gpp <- brm(log1p(GPP_ghr) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 10, iter = 5000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
fit_brs <- brm(log1p(BranchSpace_cm3) ~ DHW1  + mpa + gravtot + env1 + env2 + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 10, iter = 5000, prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
fit_ing <- brm(log1p(Accretion_kgyr) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 10, iter = 5000,  prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
fit_org <- brm(log1p(OrgGrowth_gyr) ~ DHW1  + gravtot + mpa +env1 + env2 +  (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 10, iter = 5000,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))


fit_b_env <- list(fit_cc,fit_cal, fit_sto, fit_rug, fit_gpp, fit_brs, fit_org, fit_ing)

}



#fit no env

ind_benthic_models_noenv <- function(benthos) {
  
pr <- get_prior(log1p(Calc_ghr) ~ DHW1  + gravtot + mpa + (1|SiteCode), data = benthos)
pr$prior[which(pr$class == "b")] <- "normal(0,1)"

fit_cc_noenv <- brm(log1p(coral) ~ DHW1  + gravtot + mpa + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
fit_cal_noenv <- brm(log1p(Calc_ghr) ~ DHW1  + gravtot + mpa + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
fit_sto_noenv <- brm(log1p(Storage_kg) ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000, seed = 10,  prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
fit_rug_noenv <- brm(log1p(Rugosity) ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
fit_gpp_noenv <- brm(log1p(GPP_ghr) ~ DHW1  + gravtot + mpa + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
fit_brs_noenv <- brm(log1p(BranchSpace_cm3) ~ DHW1  + mpa + gravtot  + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
fit_ing_noenv <- brm(log1p(Accretion_kgyr) ~ DHW1  + gravtot + mpa + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
fit_org_noenv <- brm(log1p(OrgGrowth_gyr) ~ DHW1  + gravtot + mpa +  (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))


fit_b_noenv <- list(fit_cc_noenv,fit_cal_noenv, fit_sto_noenv, fit_rug_noenv, fit_gpp_noenv, fit_brs_noenv, fit_org_noenv, fit_ing_noenv)

}



###########  fish models 

ind_fish_models_env <- function(fish) {

library(brms)

pr <- get_prior(log1p(herb) ~ DHW1  + gravtot + mpa + env1 + env2 +  (1|SiteCode), data = fish)
pr$prior[1] <- "normal(0,1)"


#fit env

fit_her <- brm(log1p(herb) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000, seed = 10,  prior = pr)
fit_pla <- brm(log1p(plank) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000, seed = 10,  prior = pr)
fit_pis <- brm(log1p(pisc) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000, seed = 10,  prior = pr)
fit_pro <- brm(log1p(prod) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000, seed = 10,  prior = pr)
fit_exP <- brm(log1p(exP) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000, seed = 10,  prior = pr)
fit_exN <- brm(log1p(exN) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000, seed = 10,  prior = pr)
fit_turn <- brm(log1p(turn) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000, seed = 10,  prior = pr)


fit_f_env <- list(fit_her, fit_pla, fit_pis, fit_pro, fit_exP, fit_exN, fit_turn)

}


#fit no env

ind_fish_models_noenv <- function(fish) {
  
pr <- get_prior(log1p(herb) ~ DHW1  + gravtot + mpa + (1|SiteCode), data = fish)
pr$prior[1] <- "normal(0,1)"

fit_her_noenv <- brm(log1p(herb) ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000, seed = 10,  prior = pr)
fit_pla_noenv <- brm(log1p(plank) ~ DHW1  + gravtot + mpa + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000, seed = 10,  prior = pr)
fit_pis_noenv <- brm(log1p(pisc) ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000, seed = 10,  prior = pr)
fit_pro_noenv <- brm(log1p(prod) ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000, seed = 10,  prior = pr)
fit_exP_noenv <- brm(log1p(exP) ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000, seed = 10,  prior = pr)
fit_exN_noenv <- brm(log1p(exN) ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000, seed = 10,  prior = pr)
fit_turn_noenv <- brm(log1p(turn) ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 10, iter = 1000, seed = 10,  prior = pr)


fit_f_noenv <- list(fit_her_noenv, fit_pla_noenv, fit_pis_noenv, fit_pro_noenv, fit_exP_noenv, fit_exN_noenv, fit_turn_noenv)
}





plot_pc <- function(labels2, models_env, models_noenv) {
  
  library(tidyverse)
  library(rlang)
  
  #load("output/fit_pcs.RData")
  #load("output/fit_pcs_noenv.RData")
  
  #fit_pcs <- fish_models_env
  #fit_pcs_noenv <- fish_models_noenv
  
  stat_interval <- function(..., .width = .95) {
    stat_summary(..., fun.data = ggdist::median_hdci,
                 geom = "linerange", fun.args = list(.width = .width))
  }
  stat_median <- function(...) {
    stat_summary(..., fun = median, geom = "point")
  }
  extract_named_vec <- function(df, col, names_from = "fun") {
    out <- df[[col]]
    names(out) <- df[[names_from]]
    out
  }
  ref <- labels2
  model_list <- list(
    `(1) PC1` = list(models_env[[1]], models_noenv[[1]]),
    `(2) PC2` = list(models_env[[2]], models_noenv[[2]]),
    `(3) PC3` = list(models_env[[3]], models_noenv[[3]]),
    `(4) PC4` = list(models_env[[4]], models_noenv[[4]])
  )
  sample_post <- purrr::map_dfr(model_list, function(x) {
    noenv <- x[[2]]
    x <- x[[1]]
    r2 <- brms::bayes_R2(x, summary = FALSE, re_formula = NA)
    r2_noenv <- brms::bayes_R2(noenv, summary = FALSE, re_formula = NA)
    brms::as_draws_df(x, variable = c("b_DHW1", "b_gravtot", "b_mpaNotake")) |>
      dplyr::mutate(
        `italic(R^2) * \" w/o env.\"` = .env$r2_noenv[, 1],
        `italic(R^2) * \" w/ env.\"` = .env$r2[, 1],
        b_gravtot = b_gravtot * 1e4
      ) |>
      dplyr::rename(
        DHW = b_DHW1, `\"Gravity (x 10\"^-4 * \")\"` = b_gravtot,
        `No-take` = b_mpaNotake
      ) |>
      dplyr::select(
        DHW, `\"Gravity (x 10\"^-4 * \")\"`, `No-take`,
        `italic(R^2) * \" w/ env.\"`,`italic(R^2) * \" w/o env.\"`
      ) |>
      tidyr::pivot_longer(
        cols = everything(), names_to = "predictor", values_to = "value"
      ) |>
      suppressWarnings()
  }, .id = "lab") |>
    dplyr::left_join(ref, by = join_by(lab)) |>
    dplyr::mutate(
      lab = factor(
        lab, levels = rev(c("(1) PC1", "(2) PC2", "(3) PC3", "(4) PC4"))
      ),
      predictor = factor(
        predictor, levels = c(
          "DHW", "\"Gravity (x 10\"^-4 * \")\"", "No-take",
          "italic(R^2) * \" w/o env.\"", "italic(R^2) * \" w/ env.\""
        )
      )
    )
  cols <- extract_named_vec(ref, "cols")
  labs <- extract_named_vec(ref, "lab")
  fig3_b <- ggplot(data = sample_post) +
    aes(x = lab, y = value, colour = lab) +
    geom_hline(
      data = dplyr::filter(sample_post, !grepl("env\\.", predictor)),
      mapping = aes(yintercept = 0), linetype = 2
    ) +
    stat_interval(.width = .80, linewidth = 1.6, colour = "black") +
    stat_interval(.width = .95, linewidth = 0.6, colour = "black") +
    stat_interval(.width = .80, linewidth = 1.5) +
    stat_interval(.width = .95, linewidth = 0.5) +
    stat_median(size = 1.8, colour = "black") +
    stat_median(size = 1.7) +
    scale_colour_manual(values = cols) +
    coord_flip() +
    labs(x = "Function", y = NULL) + 
    labs(caption = "Medians with 95% (thin), 80% (thick) C.I.") +
    facet_grid(~ predictor, scales = "free_x", labeller = label_parsed) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      panel.spacing = unit(0.8, "cm", data = NULL)
    )
  ggsave("output/fig3_b_pcs.png", fig3_b, dev = "png", width = 10, height = 3.2)
  
  
}








