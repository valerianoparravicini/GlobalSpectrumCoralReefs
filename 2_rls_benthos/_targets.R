library(targets)

source("R/functions.R")

dir.create("output", recursive = TRUE, showWarnings = FALSE)

options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse", "fishualize", "patchwork", "readxl", "sp", "ggforce", "scales", "reshape2", 
                            "cowplot", "ggrepel", "fishualize", "viridis", "ggcorrplot", "parallel", "utils"))



list(
  # data files
  tar_target(rls_raw_file, "data/rls_out.csv", format = "file"),
  tar_target(rls_info_file, "data/rls_info.csv", format = "file"),
  tar_target(traits_file, "data/trait_dat.csv", format = "file"),
  tar_target(geo_file, "data/rls_geo.csv", format = "file"),
  tar_target(rls_sites_file, "data/RLS_sitesInfos.rds", format = "file"),
  tar_target(rls_predictor_file, "data/RLS_predictors.csv", format = "file"),
  
  # load data
  tar_target(rls_raw, read_csv(rls_raw_file)),
  tar_target(rls_info, read_csv(rls_info_file)),
  tar_target(traits, read_csv(traits_file)),
  tar_target(geo, read_csv(geo_file)),
  tar_target(rls_sites, readRDS(rls_sites_file)),
  tar_target(pred, read_csv(rls_predictor_file)),
  
  # compute benthic functions
  tar_target(rls, 
             clean_rls_dataset(rls_raw, rls_info)),
  
  tar_target(sim_size, 
             sim_colony_size_rls(rls)),
  
  tar_target(indivs, 
             surveys_size_arranged(sim_size)),
  
  tar_target(morphs_traits, 
             geometry_change_1_year(traits, indivs, geo, rls_info)),
  
  tar_target(benthic_functions, 
             compute_functions(morphs_traits, pred, rls_info)),
  
  tar_target(taxa_contribution, 
             taxa_contribution_plot(rls_info, morphs = benthic_functions)),
  
  tar_target(data_corrected_size, 
             site_level_functions(rls, rls_sites, rls_info, morphs = benthic_functions, morphs_traits))
  
)
