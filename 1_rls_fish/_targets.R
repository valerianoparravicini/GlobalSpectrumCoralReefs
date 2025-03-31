library(targets)

source("R/functions.R")

dir.create("output", recursive = TRUE, showWarnings = FALSE)
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse", "fishualize", "patchwork", "brms"))

list(
  # data files
  tar_target(rls_fish_file, "data/RLS_fish_2021.RDS", format = "file"),
  tar_target(rls_meta_file, "data/RLS_MEOW_meta.csv", format = "file"),
  tar_target(rls_predictor_file, "data/RLS_predictors.csv", format = "file"),
  tar_target(params_file, "data/params_sst_size.rds", format = "file"),
  tar_target(troph_file, "data/extrapolation_trophic_guilds.csv", format = "file"),

  # load data
  tar_target(rls_fish, readRDS(rls_fish_file)),
  tar_target(rls_predictors, read_csv(rls_predictor_file)),
  tar_target(rls_meta, read_csv(rls_meta_file)),
  tar_target(params, readRDS(params_file)),
  tar_target(troph, read_csv(troph_file)),
  
  # select species, sizes, sst
  tar_target(sp_size_sst, 
             unique_size_sst(rls_fish, rls_predictors)),
  
  tar_target(sp_size_sst_params, 
             add_params(sp_size_sst, params)),
  
  tar_target(fluxes_sp_size_sst, 
             run_fishflux_segmented(sp_size_sst_params, cores = 4)),
  
  tar_target(fluxes_sp_size_sst_out, 
             write_csv(fluxes_sp_size_sst, "output/fluxes_sp_size_sst_rls.csv")),
  
  tar_target(functions_sp_size_sst, 
             ind_functions(fluxes_sp_size_sst, troph, sp_size_sst_params)),
  
  tar_target(transect_ind_functions, 
             combine_transects_functions(rls_fish, rls_predictors, rls_meta,  functions_sp_size_sst)),
  
  tar_target(summary_fish_functions, 
             summarize_transect(transect_ind_functions)),
  
  tar_target(summary_fish_functions_output, 
             write_csv(summary_fish_functions, "output/transect_fish_functions.csv"))
  

)
