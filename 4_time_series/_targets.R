library(targets)

source("R/functions.R")

dir.create("output", recursive = TRUE, showWarnings = FALSE)
options(tidyverse.quiet = TRUE)

tar_option_set(packages = c("tidyverse", "patchwork", "readr","tidybayes",  "brms", "fishflux", "ggplot2", "ggcorrplot"))

# data files
list(
  
  tar_target(ts_fish_file, "data/fish.time.clean.rds", format = "file"),
  tar_target(ts_benthic_file, "data/benthic.time.clean.csv", format = "file"),
  tar_target(params_file, "data/params_sst_size.rds", format = "file"),
  tar_target(troph_file, "data/extrapolation_trophic_guilds.csv", format = "file"),
  tar_target(ts_benthos_case_study_file, "data/ts_benthic_data.csv", format = "file"),
  tar_target(rls_info_file, "../2_rls_benthos/data/rls_info.csv", format = "file"),
  tar_target(traits_file, "../2_rls_benthos/data/trait_dat.csv", format = "file"),
  tar_target(geo_file, "../2_rls_benthos/data/rls_geo.csv", format = "file"),
  tar_target(predictors_ts_file, "data/predictors_ts.csv", format = "file"),
  
  
# load data
  tar_target(fish_raw, readRDS(ts_fish_file)),
  tar_target(benthic_raw, read_csv(ts_benthic_file)),
  tar_target(
    params, 
    readRDS(params_file) %>%
      group_by(Family, Genus, Species, species, diet, v_m) %>%
      summarize_all(mean) %>%
      ungroup()
  ),
  tar_target(troph, read_csv(troph_file)),
  tar_target(fish, clean_fish(fish_raw)),
  tar_target(coral, clean_bent(benthic_raw)),
  tar_target(sst_list,sst_time_series(fish)),
  tar_target(sst_data, 
    lapply(sst_list, function(x){attributes(x)$info %>% mutate(sst = x$sst)}) %>% plyr::ldply()),
  tar_target(ts_benthos_case_study, read_csv(ts_benthos_case_study_file)),
  tar_target(rls_info, read_csv(rls_info_file)),
  tar_target(traits, read_csv(traits_file)),
  tar_target(geo, read_csv(geo_file)),
  tar_target(predictors_ts, read_csv(predictors_ts_file)),
  
  # select species, sizes, sst
  tar_target(sp_size_sst, unique_spsize_sst(fish, sst_data)),
  tar_target(sp_size_params, add_params(sp_size_sst, params)),
  tar_target(fluxes_sp_size, run_fishflux_segmented(sp_size_params, cores = 10)),
  tar_target(functions_sp_size, ind_functions(fluxes_sp_size, troph, sp_size_params)),
  tar_target(transect_ind_functions,combine_transects_functions(fish, functions_sp_size)),
  tar_target(summary_fish_functions, summarize_transect(transect_ind_functions, coral)),
  tar_target(summary_fish_functions_output, write_csv(summary_fish_functions, "output/ts_transect_fish_functions.csv")),
  tar_target(ts_benthos_case_study_clean, clean_benthos_ts(ts_benthos_case_study)),
  tar_target(indivs, sim_colony_size_ts(ts_benthos_case_study = ts_benthos_case_study_clean)),
  tar_target(morphs_traits, geometry_change_1_year_ts(traits, indivs, geo, rls_info)),
  tar_target(benthic_functions, compute_functions(morphs_traits, predictors_ts, rls_info, ts_benthos_case_study = ts_benthos_case_study_clean)),
  tar_target(ts_functions_dataset, create_ts_functions_dataset(benthos = benthic_functions, fish = summary_fish_functions)),
  tar_target(ts_pca, predict_pca(ts_functions_dataset)),
  tar_target(Fig_4, fig_4(ts_pca))
)
  

