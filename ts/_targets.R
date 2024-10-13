library(targets)
source("R/functions.R")
dir.create("output", recursive = TRUE, showWarnings = FALSE)
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse", "patchwork", "readr","tidybayes",  "brms", "fishflux"))
list(
  # data files
  tar_target(
    ts_fish_file, "data/fish.time.clean.rds", format = "file"
  ),
  tar_target(
    ts_benthic_file, "data/benthic.time.clean.csv", format = "file"
  ),
  tar_target(
    params_file, "../rls/data/params_sst_size.rds", format = "file"
  ),
  tar_target(
    troph_file, "../rls/data/extrapolation_trophic_guilds.csv", format = "file"
  ),
  # load data
  tar_target(
    fish_raw, readRDS(ts_fish_file)
    ),
  tar_target(
    benthic_raw, read_csv(ts_benthic_file)
  ),
  tar_target(
    params, 
    readRDS(params_file) %>%
      group_by(Family, Genus, Species, species, diet, v_m) %>%
      summarize_all(mean) %>%
      ungroup()
  ),
  tar_target(
    troph, read_csv(troph_file)
  ),
  tar_target(
    fish, clean_fish(fish_raw)
  ),
  tar_target(
    coral, clean_bent(benthic_raw)
  ),
  tar_target(
    sst_list,
    sst_time_series(fish)),
  tar_target(
    sst_data, 
    lapply(sst_list, function(x){attributes(x)$info %>% mutate(sst = x$sst)}) %>% plyr::ldply()),

  # select species, sizes, sst
  tar_target(
    sp_size_sst, unique_spsize_sst(fish, sst_data)
  ),
  tar_target(
    sp_size_params, add_params(sp_size_sst, params)
  ),
  tar_target(
    fluxes_sp_size, run_fishflux_segmented(sp_size_params, cores = 5)
  ),
  tar_target(
    functions_sp_size, ind_functions(fluxes_sp_size, troph, sp_size_params)
  ),
  tar_target(transect_ind_functions,
             combine_transects_functions(fish, functions_sp_size)),
  tar_target(summary_fish_functions, summarize_transect(transect_ind_functions, coral)),
  tar_target(summary_fish_functions_output, write_csv(summary_fish_functions, "output/ts_transect_fish_functions.csv"))
  
)
