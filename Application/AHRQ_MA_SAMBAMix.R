library(tidyverse)
library(nimble)
source("~/bayesmbmm/R/helpers.R")
source("~/bayesmbmm/R/bayesmbmm_OFMM_SF_v2.R")

# 
# args <- commandArgs(trailingOnly = TRUE)
# job_id <- as.numeric(args[1])

set.seed(2026)
# #Read in AHRQ data
ahrq<-readRDS(file = "~/AHRQData/sdoh_ahrq_statesCTs_complete.rds")
#names(ahrq)<- c("NY", "MA")

madat<-ahrq[[2]] 
#summary(madat)

start_time <- Sys.time()
# fit<-Ombmm_FeatSaliency(
#   x = madat,
#   K = 30,
#   niter = 20000,
#   nburnin = 7000,
#   nchain = 1,
#   init_list = NULL,
#   g =  2, h = 5, 
#   a_mu0 = 2, b_mu0 = 8, #more vars near 0, so this prior anchors irrelevant features near realistic values
#   a_mu1 = 2, b_mu1 = 2,
#   a_rho = 1, b_rho = 8, #mild sparsity; variables near 0 will def be noisy
#   e0 = 0.1,
#   priorOnAlpha = "no", plot_trace = FALSE)

fit <- Ombmm_FeatSaliency(
  x = madat,
  K = 30,
  niter = 20000,
  nburnin = 10000,
  nchain = 1,
  # Increase Cluster Density
  g = 15, h = 1,        # Mean=15, tighter clusters help separate features
  a_mu0 = 1, b_mu0 = 1, # Flat prior for irrelevant features. 
  # This prevents the model from 'dumping' variables 
  # into mu0 just because they are near zero.
  #  Increase Feature Inclusion Pressure
  a_rho = 2, b_rho = 5, # Encourages more features (Expected inclusion ~28%)
  #Force Component Sparsity
  e0 = 0.001,           # results: ahrqMA_sambo_3.27.26.rds
  #e0 = 0.01,       # results:ahrqMA_sambo_3.28.26.rds
  priorOnAlpha = "no", 
  plot_trace = FALSE
)

end_time <- Sys.time()
runtime_sec <- as.numeric(difftime(end_time, start_time, units = "secs")) 

#out_file <- sprintf("/n/netscratch/stephenson_lab/Lab/crodriguez/AHRQ_application/ahrqMA_sambamix_8.24.26.rds")
out_file <- sprintf("/n/netscratch/stephenson_lab/Lab/crodriguez/AHRQ_application/ahrqMA_sambamix_8.28.26.rds")

saveRDS(list(fit, runtime_sec), out_file)





