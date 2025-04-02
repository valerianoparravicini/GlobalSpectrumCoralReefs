# Code for "The global spectrum of coral reef ecosystem functions in the Anthropocene"

This repository contains a number of functions to reproduce the analysis presented in "The global spectrum of coral reef ecosystem functions in the Anthropocene".
The repository contains 4 distinct folders, and each folder cotains a separate code pipeline organized using the 'target' package in R.
There are 4 pipelines:

1) 1_rls_fish -> Computes fish functions on the Reef Life Survey dataset;
2) 2_rls_benthos -> Computes benthic functions on the Reef Life Survey dataset;
3) 3_rls_analysis -> Performs spatial analyses, bayesian models and produces Fig1, Fig2 and Fig3;
4) 4_time_series -> Computes fish and benthic functions on the time series datasets and produces Fig. 4.

To reproduce the analyses it is necessary to run targets::tar_make (i.e. run the pipeline) for each repository in the order they have been presented above.

# NOTES: 
1) Each pipeline has an associates renv environment produced with the 'renv' package in R to ensure the installation of the exact version of R packages used to conduct the analyses. 
2) The analyses have been performed with the "R version 4.4.2 (2024-10-31)" installed on a Debian GNU/Linux 12 (bookworm) machine.
:exclamation: :boom: **CAUTION**: The script depends on parallel
computation and uses up to **120 threads** and **250G** of memory, and thus
should be run on a supercomputer. It takes about one four days to reproduce
the entire project.

# sessionInfo()

R version 4.4.2 (2024-10-31)
Platform: x86_64-pc-linux-gnu
Running under: Debian GNU/Linux 12 (bookworm)

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.11.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.11.0

locale:
 [1] LC_CTYPE=fr_FR.UTF-8       LC_NUMERIC=C               LC_TIME=fr_FR.UTF-8        LC_COLLATE=fr_FR.UTF-8     LC_MONETARY=fr_FR.UTF-8   
 [6] LC_MESSAGES=fr_FR.UTF-8    LC_PAPER=fr_FR.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/Paris
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.4.2    tools_4.4.2       rstudioapi_0.17.1 renv_1.1.4   
