# Code for "The global spectrum of coral reef ecosystem functions in the Anthropocene"

This repository contains a number of functions to reproduce the analysis presented in "The global spectrum of coral reef ecosystem functions in the Anthropocene".
The repository contains 4 distinct folders, and each folder cotains a separate code pipeline organized using the 'target' package in R.
These are 4 pipelines:

1) rls -> quantify fish-associated functions on the RLS dataset
2) ts -> quantify fish-associated functions on the time series dataset
3) multi -> quantify benthic functions and performs multivariate analyses (the code produces also Fig1, Fig2, Fig3ab, Fig4)
4) models -> perform the bayesian models and produces Fig3c

To reproduce the analyses it is necessary to run targets::tar_make (i.e. run the pipeline) for each repository in the order they have been presented above.

NOTES: 
1) the models pipeline is coded to run in parallel on 40 threads and takes ~2 days to run, or more, according to the processor
2) the code requires the following R packages : library(readxl) ; library(sp) ; library(ggforce) ; library(tidyverse); library(patchwork) ; library(scales)
  library(reshape2) ; library(cowplot) ; library(ggrepel) ; library(fishualize) ; library(viridis) ; library(ggcorrplot); library(sf) ; library(brms) ; library(fishflux)  ; library(rlang) 
