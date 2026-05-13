#=================================================================
# Use simulated data from ACSmock (see R/simulate_data_cases2.R)
#For each case,  we generate 100 datasets and calculate MSE for all model parameters
# Programmer: Carmen Rodriguez C.
# Date updated: 2/6/2026
#===========================================================

# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
job_id<-as.numeric(args[1])

# Load libraries
setwd("/n/home03/crodriguezcabrera/bayesmbmm/")
library(nimble)
library(coda)
library(tictoc)
library(tidyverse)
library(dplyr)
source("R/helpers.R")
source("R/bayesmbmm_OFMM_SF_v2.R")

set.seed(2025 + job_id)

datlist7<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_case7.rds") 
datlist5<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_case5.rds")


alloclist7<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case7.rds")
alloclist5<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case5.rds")
alloclist6<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case6.rds")

sim_data<- list(case7 = datlist7, case5= datlist5)
#list of true allocations per scenario
alloc_true_all<- list(case7=alloclist7, case5=alloclist5)
scenariosdat<- alloclist6[[1]] #tibble of true parameters (mu, lambda, pi)
scenariosdat<-scenariosdat %>% rename(dgp_id = cid)


# Define all scenario × replicate pairs
scenarios <- names(sim_data)  
jobs <- expand.grid(scen = scenarios, rep = 1:100, stringsAsFactors = FALSE)
cid <- jobs$scen[job_id]
rid <- jobs$rep[job_id]


alloc_true <- alloc_true_all[[cid]][[rid]]
true_params<- get_true_params(scenariosdat, cid)
p<-true_params$numfeat

casedat<-sim_data[[cid]]

#Save result 
out_file <- sprintf("/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/Ombmmcases_old_%s_rep%d.rds", cid, rid)

start_time <- Sys.time()
# Fit model
fit <- tryCatch(
  {Ombmm_FeatSaliency(
      x = casedat[[rid]],
      K = 20,
      niter = 15000, 
      nburnin = 5000,
      nchain = 1,
      init_list = NULL, 
      g =  2, h = 1, 
      a_mu0 = 2, b_mu0 = 6, 
      a_mu1 = 2, b_mu1 = 2,
      a_rho = 1, b_rho = 5, 
      e0=0.1, priorOnAlpha = "no",
      p.true = true_params$true_pi,
      plot_trace = FALSE)
    
  },
  error = function(e) {
    message(sprintf("Job %d | Dataset %d failed: %s", job_id, rid, e$message))
    return(NULL)
  }
)
end_time <- Sys.time()
runtime_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))

if (is.null(fit)) {
  saveRDS(
    list(
      job = job_id,
      cid = cid,
      Kmode = NA,
      mse_lambda = NA,
      mab_lambda = NA,
      MSE_mu = NA,
      ARI = NA,
      runtime = runtime_sec,
      failed = TRUE
    ),
    out_file
  )
  quit(save = "no", status = 0)
}


#Posterior summaries
sSum <- samplesSummary(fit$modeloutput_perchain$parameters.ECR.mcmc)
Kmode<- fit$K
alloc_est<-fit$modeloutput_perchain$clusterMembershipPerMethod[1,]
lambda_est <- as.numeric(sSum[paste0("lambda[", 1:Kmode, "]"), 1])
if (length(lambda_est) != length(true_params$true_lambda)) {
  mse_lambda <- NA
} else {
  mse_lambda<-get_dist(lambda_est, true_params$true_lambda, dist_type = "mean_sq")
}

mu_mat_est <- sapply(1:p, function(j) {
  as.numeric(sSum[paste0("mu[", 1:Kmode, ", ", j, "]"), 1])})

#Check on mu and ARI
true_mu<- true_params$true_mu
if (nrow(mu_mat_est) != nrow(true_mu) ||
    ncol(mu_mat_est) != ncol(true_mu)) {
  MSE_mu<- NA
} else {
  MSE_mu <- get_dist(mu_mat_est, true_mu, dist_type = "mean_sq")
}

#MSE_mu_vec[d]<-get_dist(mu_mat_est, true_mu, dist_type = "mean_sq")

# ARI
if (length(alloc_est) != length(alloc_true)) {
  ARI <- NA
} else {
  ARI <- mclust::adjustedRandIndex(alloc_est, alloc_true)
}



saveRDS(
  list(
    job_id = job_id,
    cid = cid,
    rid = rid,
    runtime = runtime_sec,
    Kmode = Kmode,
    lambda_est = lambda_est, 
    mse_lambda = mse_lambda, 
    mu_est = mu_mat_est, 
    MSE_mu = MSE_mu,
    ARI = ARI,
    nk = fit$nk_iter,
    samplesall = fit$samples_raw,
    PIP = as.vector(fit$FeatureSal_out[[2]])),
  out_file
)





















# 
# res<- run_1sim_OmbmmSF(cid, rid, sim_data, true_alloc, scenariosdat,
#                        Kmax = 20, niter = 15000, nburnin = 5000,
#                        g =  2, h = 1, a_mu0 = 2, b_mu0 = 6, a_mu1 = 2, b_mu1 = 2,
#                        a_rho = 1, b_rho = 5, e0=0.1, priorOnAlpha = "no")

# 
# end_time <- Sys.time()
# runtime_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))
# res$runtime_sec <- runtime_sec
# 
# 
# #Save result 
# out_file <- sprintf("/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/OmbmmresultACSmock_cniter1_%s_rep%d.rds", cid, rid)
# saveRDS(res, out_file)
# #Submitted batch job 64332011 (e0 = 0.01)
# 
# 
# #3/21/26
# 
# true_pi<-scenariosdat$pivec[[1]]
# 
# casedat<-datlist5[[25]]
# fit <-Ombmm_FeatSaliency(
#       x = casedat,
#       K = 5,
#       niter = 5000, 
#       nburnin = 2000,
#       nchain = 1,
#       init_list = NULL, 
#       g =  5, h = 3, 
#       a_mu0 = 1, b_mu0 = 1, 
#       a_mu1 = 1, b_mu1 = 1,
#       a_rho = 1, b_rho = 1, 
#       e0=10, priorOnAlpha = "no",
#       p.true = true_pi,
#       plot_trace = TRUE
#     )
# 
# 
# 
# 
# 
# 
# 
# 
