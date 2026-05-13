library(tidyverse)
library(nimble)
source("R/helpers.R")
source("R/bayesmbmm_OFMM_SF_v2.R")
# 
args <- commandArgs(trailingOnly = TRUE)
job_id <- as.numeric(args[1])
# 
# hp_all <- readRDS("~/AHRQData/hp_grid.rds")
# 
# this_hp <- hp_all[job_id, ]
# e0 <-this_hp$e0
# a_rho<-this_hp$a_rho
# b_rho<-this_hp$b_rho
# sid<- this_hp$dataset
# 
# 
# set.seed(1000 + job_id)
# 
# #Read in AHRQ data
# ahrq<-readRDS(file = "~/AHRQData/sdoh_ahrq_statesCTs_complete.rds")
# names(ahrq)<- c("NY", "MA")
# 
# 
# start_time <- Sys.time()
# 
# fit<-tryCatch({
#   Ombmm_FeatSaliency(
#   x = ahrq[[sid]],
#   K = 20,
#   niter = 15000,
#   nburnin = 5000,
#   nchain = 1,
#   init_list = NULL, 
#   g =  5, h = 3, a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
#   a_rho = a_rho, b_rho = b_rho, e0 = e0,
#   priorOnAlpha = "no", plot_trace = FALSE)
# }, error = function(e) {
#   message("failed: ", e$message)
#   return(NULL)
# })
# 
# end_time <- Sys.time()
# runtime_sec <- as.numeric(difftime(end_time, start_time, units = "secs")) 
# # 
# outfile <- paste0(
#   "AHRQ_application/job_",
#   job_id,
#   "_e0_", this_hp$e0,
#   "_rho_", this_hp$rho_mean,
#   "_dat_", this_hp$dataset,
#   ".rds"
# )
# 
# saveRDS(
#   list(
#     job_id = job_id,
#     hyperparameters = this_hp,
#     runtime_sec = runtime_sec,
#     fit= fit
#   ),
#   outfile
# )

set.seed(2026 + job_id)
# #Read in AHRQ data
ahrq<-readRDS(file = "~/AHRQData/sdoh_ahrq_statesCTs_complete.rds")
names(ahrq)<- c("NY", "MA")

sid<- names(ahrq)[job_id]

start_time <- Sys.time()
fit<-Ombmm_FeatSaliency(
  x = ahrq[[sid]],
  K = 30,
  niter = 5000,
  nburnin = 0,
  nchain = 1,
  init_list = NULL,
  g =  5, h = 3, a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
  a_rho = 1, b_rho = 1, e0 = 0.01,
  priorOnAlpha = "no", plot_trace = FALSE)


end_time <- Sys.time()
runtime_sec <- as.numeric(difftime(end_time, start_time, units = "secs")) 

out_file <- sprintf("/n/netscratch/stephenson_lab/Lab/crodriguez/AHRQ_application/ahrq_%s.rds", sid)

saveRDS(list(fit, runtime_sec), out_file)

