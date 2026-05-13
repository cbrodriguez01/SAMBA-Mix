library(tidyverse)
library(nimble)
source("~/bayesmbmm/R/helpers.R")
source("~/bayesmbmm/R/bayesmbmm_OFMM_SF_v2.R")

set.seed(2026)
# #Read in AHRQ data
ahrq<-readRDS(file = "~/AHRQData/sdoh_ahrq_statesCTs_complete.rds")
nydat<-ahrq[[1]]

start_time <- Sys.time()
fit<-Ombmm_FeatSaliency(
  x = nydat,
  K = 30,
  niter = 20000,
  nburnin = 10000,
  nchain = 1,
  init_list = NULL,
  # Increase Cluster Density
  g = 2, h = 1,        # Mean=15, tighter clusters help separate features
  a_mu0 = 1, b_mu0 = 1, # Flat prior for irrelevant features. 
  # This prevents the model from 'dumping' variables 
  # into mu0 just because they are near zero.
  #  Increase Feature Inclusion Pressure
  a_rho = 1, b_rho = 15, # Encourages more features (Expected inclusion ~28%)
  #Force Component Sparsity
  #e0 = 0.001,          
  e0 = 0.001,       
  priorOnAlpha = "no", 
  plot_trace = FALSE)

end_time <- Sys.time()
runtime_sec <- as.numeric(difftime(end_time, start_time, units = "secs")) 

out_file <- sprintf("/n/netscratch/stephenson_lab/Lab/crodriguez/AHRQ_application/ahrqNY_sambo_4.8.26.rds")

saveRDS(list(fit, runtime_sec), out_file)























# /ahrqNY_sambo_4.4.26.rds
# fit<-Ombmm_FeatSaliency(
#   x = nydat,
#   K = 30,
#   niter = 20000,
#   nburnin = 10000,
#   nchain = 1,
#   init_list = NULL,
#   # Increase Cluster Density
#   g = 15, h = 1,        # Mean=15, tighter clusters help separate features
#   a_mu0 = 1, b_mu0 = 1, # Flat prior for irrelevant features. 
#   # This prevents the model from 'dumping' variables 
#   # into mu0 just because they are near zero.
#   #  Increase Feature Inclusion Pressure
#   a_rho = 2, b_rho = 5, # Encourages more features (Expected inclusion ~28%)
#   #Force Component Sparsity
#   #e0 = 0.001,          
#   e0 = 0.0002,       
#   priorOnAlpha = "no", 
#   plot_trace = FALSE)
# 
# “Make clusters tight”
# “Use very few clusters”
# “Don’t strongly penalize feature inclusion”
# 



# delta <- abs(mu_matrix[1,] - mu_matrix[2,])
# sort(delta, decreasing = TRUE)


# #I could do this together but i know MA is going to finish first so i'll separate
# 
# library(tidyverse)
# library(nimble)
# source("~/bayesmbmm/R/helpers.R")
# source("~/bayesmbmm/R/bayesmbmm_OFMM_SF_v2.R")
# 
# args <- commandArgs(trailingOnly = TRUE)
# job_id <- as.numeric(args[1])
# 
# set.seed(2026 + job_id)
# # #Read in AHRQ data
# ahrq<-readRDS(file = "~/AHRQData/sdoh_ahrq_statesCTs_complete.rds")
# names(ahrq)<- c("NY", "MA")
# 
# hp_all <- readRDS("~/AHRQData/hp_grid.rds")
# 
# this_hp <- hp_all[job_id, ]
# e0 <-this_hp$e0
# a_rho<-this_hp$a_rho
# b_rho<-this_hp$b_rho
# 
# start_time <- Sys.time()
# fit<-Ombmm_FeatSaliency(
#   x = ahrq[["NY"]],
#   K = 25,
#   niter = 4000,
#   nburnin = 2000,
#   nchain = 1,
#   init_list = NULL,
#   g =  5, h = 3, a_mu0 = 1, b_mu0 = 1,a_mu1 = 1, b_mu1 = 1,
#   a_rho = a_rho, b_rho = b_rho, e0 = e0,
#   priorOnAlpha = "no", plot_trace = FALSE)
# 
# 
# end_time <- Sys.time()
# runtime_sec <- as.numeric(difftime(end_time, start_time, units = "secs")) 
# 
# out_file <- sprintf("/n/netscratch/stephenson_lab/Lab/crodriguez/AHRQ_application/ahrqNY_model_%s.rds", job_id)
