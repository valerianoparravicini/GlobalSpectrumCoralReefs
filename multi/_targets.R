library(targets)
source("R/functions.R")

dir.create("output", recursive = TRUE, showWarnings = FALSE)

options(tidyverse.quiet = TRUE)

tar_option_set(packages = c("tidyverse", "brms", "readxl"))


list(
  # data files
  tar_target(data_ts_coord_file, "data/data_ts_coord.xlsx", format = "file"),
  tar_target(coord_case_study_file, "data/coord_ts_case_study.xlsx", format = "file"), 
  tar_target(rls_file, "data/rls_out.csv", format = "file"),
  tar_target(rls_info_file, "data/rls_info.csv", format = "file"),
  tar_target(pred_file, "data/RLS_predictors.csv", format = "file"),
  tar_target(rls_sites_file, "data/RLS_sitesInfos.rds", format = "file"),
  tar_target(geo_file, "data/rls_geo.csv", format = "file"),
  tar_target(traits_file, "data/trait_dat.csv", format = "file"), 
  tar_target(predictors_ts_file, "data/predictors_ts.csv", format = "file"),
  tar_target(ts_benthos_file, "data/03-merge_benthos.xlsx", format = "file"),
  tar_target(Seychelles_data_file, "data/Seychelles_benthic_line intercept.xlsx", format = "file"), 
  tar_target(Moorea_data_file, "data/MPA.xlsx", format = "file"),
  tar_target(ts_fish_file, "data/03-merge_fish.rds", format = "file"),
  tar_target(Benthos_Functions_file, "data/RLS_benthic_functions.csv", format = "file"),
  tar_target(Fish_Functions_file, "../rls/output/transect_fish_functions.csv", format = "file"),
  tar_target(case_study_file, "data/final_dataset.xlsx", format = "file"),
  tar_target(indivs_file, "data/rls_individuals.csv", format = "file"), 
  tar_target(Fish_functions_file, "../ts/output/ts_transect_fish_functions.csv", format = "file"),
  
  # load data
  tar_target(data_ts_coord, read_excel(data_ts_coord_file)),
  tar_target(coord_case_study, read_excel(coord_case_study_file)), 
  tar_target(rls, read_csv(rls_file)),
  tar_target(rls_info, read_csv(rls_info_file)),
  tar_target(pred, read_csv(pred_file)),
  tar_target(rls_sites, readRDS(rls_sites_file)),
  tar_target(geo, read_csv(geo_file)),
  tar_target(traits, read_csv(traits_file)), 
  tar_target(predictors_ts, read_csv(predictors_ts_file)),
  tar_target(ts_benthos, read_excel(ts_benthos_file)),
  tar_target(Seychelles_data, read_excel(Seychelles_data_file)), 
  tar_target(Moorea_data, read_excel(Moorea_data_file)),
  tar_target(ts_fish, readRDS(ts_fish_file)),
  tar_target(Benthos_Functions, read_csv(Benthos_Functions_file)),
  tar_target(Fish_Functions, read_csv(Fish_Functions_file)),
  tar_target(case_study, read_excel(case_study_file)),
  tar_target(indivs, read_csv(indivs_file)),
  tar_target(Fish_functions, read_csv(Fish_functions_file)),
  
  #compute benthic functiona
  tar_target(benthic_functions, benthic_functions_estimation(rls, traits, geo, pred, rls_info, rls_sites, indivs)),
  
  #plot
  tar_target(figure_2, plot_figure2(Fish_Functions, Benthos_Functions, case_study)),
  tar_target(figure_1, plot_figure1(Benthos_Functions,Fish_Functions, rls, data_ts_coord)),
  tar_target(figure_4, time_series(ts_benthos, ts_fish, Seychelles_data, Moorea_data, rls_info, traits, geo, predictors_ts, Fish_functions, case_study,
                                   Benthos_Functions,Fish_Functions, pred))
  
  #tar_target(dat_fish, prepare_fish_models(pred, fish)),
  
  #### models PC
  #tar_target(models_env, run_models_env(dat_mod)),
  #tar_target(models_noenv, run_models_noenv(dat_mod)),
  
  #### models ind
  #tar_target(benthic_models_env, ind_benthic_models_env(dat_benthos)),
  #tar_target(benthic_models_noenv, ind_benthic_models_env(dat_benthos)),
  #tar_target(fish_models_env, ind_fish_models_env(dat_fish)),
  #tar_target(fish_models_noenv, ind_fish_models_env(dat_fish))
  
  
  ##### plots
  #  tar_target(plot1, make_pca(transect_fish_benthic, rls_predictors))
  
)


