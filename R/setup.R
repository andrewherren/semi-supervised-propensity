#=============================================
# R setup for this project
#=============================================

## Clear .Rdata
# rm(list=ls())

## Packages
# List of packages to be installed if necessary
pkg.list <- c("bcf", "grf", "ranger", "survey", "stats", "dplyr", "magrittr",
              "assertthat", "testthat", "quantreg", "MASS", "CausalGAM", 
              "KernSmooth", "sampling", "lavaan.survey", "tmle", "sandwich", 
              "Matching", "glmnet", "dbarts", "XBART", "xtable", "here", 
              "foreach", "doParallel")
# install.packages(pkg.list)
lapply(pkg.list, require, character.only = TRUE)

## Options / configuration
# ....
# get_XBART_params <- function(y) {
#   XBART_params = list(num_trees = 15, # number of trees 
#                       num_sweeps = 60, # number of sweeps 
#                       # (samples of the forest)
#                       n_min = 1, # minimal node size
#                       alpha = 0.95, # BART prior parameter 
#                       beta = 1.25, # BART prior parameter
#                       mtry = 3, # number of variables sampled in each split
#                       burnin = 15,
#                       no_split_penality = "Auto"
#   ) # burnin of MCMC sample
#   num_tress = XBART_params$num_trees
#   XBART_params$max_depth = 250
#   XBART_params$num_cutpoints = 50;
#   # number of adaptive cutpoints
#   # prior variance of mu (leaf parameter)
#   XBART_params$tau = var(y) / num_tress 
#   return(XBART_params)
# }
# parl = FALSE
# verbose = FALSE
# dcat = 0

