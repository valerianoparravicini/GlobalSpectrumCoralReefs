# Code for "The global spectrum of coral reef ecosystem functions in the Anthropocene"

This repository contains a number of functions to reproduce the analysis presented in "The global spectrum of coral reef ecosystem functions in the Anthropocene".
The repository contains 4 distinct folders, and each folder cotains a separate code pipeline organized using the 'target' package in R.
There are 4 pipelines:

1) 1_rls_fish -> Computes fish functions on the Reef Life Survey dataset;
2) 2_rls_benthos -> Computes benthic functions on the Reef Life Survey dataset;
3) 3_rls_analysis -> Performs spatial analyses, bayesian models and produces Fig1, Fig2 and Fig3;
4) 4_time_series -> Computes fish and benthic functions on the time series datasets and produces Fig. 4.

To reproduce the analyses it is necessary to run targets::tar_make (i.e. run the pipeline) for each repository in the order they have been presented above.

NOTES: 
1) Each pipeline has an associates renv environment produced with the 'renv' package in R to ensure the installation of the exact version of R packages used to conduct the analyses. 
2) The analyses have been performed with the "R version 4.4.2 (2024-10-31)" installed on a Debian GNU/Linux 12 (bookworm) machine.
3) The 3_rls_analysis code runs using 4 cores and 30 threads per core
