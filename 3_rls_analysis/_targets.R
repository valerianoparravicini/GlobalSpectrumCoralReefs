library(targets)

source("R/functions.R")

dir.create("output", recursive = TRUE, showWarnings = FALSE)

options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse", "fishualize", "patchwork", "readxl", "sp", "ggforce", "scales", "reshape2", 
                            "cowplot", "ggrepel", "fishualize", "viridis", "ggcorrplot", "parallel", "utils", "brms", "cmdstanr",
                            "rlang", "RColorBrewer", "sf"))



list(
  # data files
  tar_target(benthic_functions_file, "../2_rls_benthos/output/RLS_benthic_functions.csv", format = "file"),
  tar_target(fish_functions_file, "../1_rls_fish/output/transect_fish_functions.csv", format = "file"),
  tar_target(rls_predictor_file, "data/RLS_predictors.csv", format = "file"),
  tar_target(labels2_file, "data/labels2.csv", format = "file"),

  

  # load data
  tar_target(rls_benthic_functions, read_csv(benthic_functions_file)),
  tar_target(fish_functions, read_csv(fish_functions_file)),
  tar_target(pred, read_csv(rls_predictor_file)),
  tar_target(labels2, read_csv(labels2_file)),


  # spatial analysis
  tar_target(list_dataset_functions, 
             combine_benthic_fish_functions(rls_benthic_functions, fish_functions)),
  
  tar_target(fig_1, 
             figure_1(list_dataset_functions)),
  
  tar_target(list_pca_functions, 
             pca_functions(list_dataset_functions, benth=rls_benthic_functions, fish=fish_functions)),
  
  tar_target(fig_2, 
             figure_2(list_pca_functions, benth=rls_benthic_functions, fish=fish_functions)),
  
  tar_target(pca, save(list_pca_functions, file = "output/pca_functions.RData")),
  
  tar_target(data_models_PC, 
             prepare_PC_data(list_pca_functions, pred)),
  
  tar_target(data_models_benthos, prepare_benthic_models(pred, benthos = rls_benthic_functions)),
  
  tar_target(data_models_fish, prepare_fish_models(pred, fish = fish_functions)),
  
  tar_target(models_env, run_models_env(dat_mod = data_models_PC)),
  
  tar_target(models_noenv, run_models_noenv(dat_mod = data_models_PC)),
  
  tar_target(benthic_models_env, ind_benthic_models_env(benthos = data_models_benthos)),
  
  tar_target(benthic_models_noenv, ind_benthic_models_noenv(benthos = data_models_benthos)),
  
  tar_target(fish_models_env, ind_fish_models_env(fish = data_models_fish)),
  
  tar_target(fish_models_noenv, ind_fish_models_noenv(fish = data_models_fish)),
  
  tar_target(fig_3c, plot_pc(labels2, models_env, models_noenv )),
  
  tar_target(fig_3ab, figure_3(list_pca_functions, pred)),
  
  tar_target(fig_3, combine_fig_3(fig_3ab, fig_3c))
)
