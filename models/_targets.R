library(targets)
source("R/functions.R")

dir.create("output", recursive = TRUE, showWarnings = FALSE)

options(tidyverse.quiet = TRUE)

tar_option_set(packages = c("tidyverse", "brms"))

list(
  # data files
  tar_target(rls_predictor_file, "../rls/data/RLS_predictors.csv", format = "file"),
  tar_target(dat_file, "data/pca_df_RLS.csv", format = "file"), 
  tar_target(benthos_file, "../multi/data/RLS_benthic_functions.csv", format = "file"),
  tar_target(fish_file, "../rls/output/transect_fish_functions.csv", format = "file"),
  tar_target(labels2_file, "data/labels2.csv", format = "file"),
  
  # load data
  tar_target(pred, read_csv(rls_predictor_file)),
  tar_target(dat, read_csv(dat_file)),
  tar_target(benthos, read_csv(benthos_file)),
  tar_target(fish, read_csv(fish_file)),
  tar_target(labels2, read_csv(labels2_file)),
  
  #prepare data PC
  tar_target(dat_mod, prepare_PC_data(dat, pred)),
  
  #prepare data ind
  tar_target(dat_benthos, prepare_benthic_models(pred, benthos)),
  tar_target(dat_fish, prepare_fish_models(pred, fish)),
  
  #### models PC
  tar_target(models_env, run_models_env(dat_mod)),
  tar_target(models_noenv, run_models_noenv(dat_mod)),
  
  #### models ind
  tar_target(benthic_models_env, ind_benthic_models_env(dat_benthos)),
  tar_target(benthic_models_noenv, ind_benthic_models_env(dat_benthos)),
  tar_target(fish_models_env, ind_fish_models_env(dat_fish)),
  tar_target(fish_models_noenv, ind_fish_models_env(dat_fish)),
  
  
  ##### plots
  tar_target(Figure_3c, plot_pc(labels2, models_env, models_noenv ))
  
)
