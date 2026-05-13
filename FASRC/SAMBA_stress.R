#=================================================================
# Use simulated data from ~/bayesmbmm/data/hardsim.R
# 50 Reps
# Programmer: Carmen Rodriguez C.
# Date updated: 3.18.26
#===========================================================

# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
rid<-as.numeric(args[1])

# Load libraries
setwd("/n/home03/crodriguezcabrera/bayesmbmm/")
library(nimble)
library(coda)
library(tictoc)
library(tidyverse)
library(dplyr)
source("R/helpers.R")
source("R/bayesmbmm_OFMM_SF_v2.R")
#out<-list(datlist, alloclist, K, mu, lambda_vec, pitrue)


output<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistStress.rds") 
datlist<-output[[1]]
alloc_true <- output[[2]]
  

p<-20
true_lambda<-output[[5]]
true_mu<-output[[4]]
true_pi<-output[[6]]


#Save result 
out_file <- sprintf("/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/SAMBO_stress2_rep%d.rds", rid)

casedat<-datlist[[rid]]

start_time <- Sys.time()
# Fit model
fit <- tryCatch(
  {
    Ombmm_FeatSaliency(
      x = casedat,
      K = 15,
      niter = 15000, 
      nburnin = 7000,
      nchain = 1,
      init_list = NULL, 
      g =  2, h = 1, 
      a_mu0 = 2, b_mu0 = 6, 
      a_mu1 = 2, b_mu1 = 2,
      a_rho = 1, b_rho = 5, 
      e0=0.1, priorOnAlpha = "no",
      p.true = true_pi,
      plot_trace = FALSE
    )
  },
  error = function(e) {
    message(sprintf("Dataset %d failed: %s", rid, e$message))
    return(NULL)
  }
)
end_time <- Sys.time()
runtime_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))

if (is.null(fit)) {
  saveRDS(
    list(
      rid = rid,
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
if (length(lambda_est) != length(true_lambda)) {
  mse_lambda <- NA
} else {
  mse_lambda<-get_dist(lambda_est, true_lambda, dist_type = "mean_sq")
}

mu_mat_est <- sapply(1:p, function(j) {
  as.numeric(sSum[paste0("mu[", 1:Kmode, ", ", j, "]"), 1])})

#Check on mu and ARI
#true_mu<- true_params$true_mu
if (nrow(mu_mat_est) != nrow(true_mu) ||
    ncol(mu_mat_est) != ncol(true_mu)) {
  MSE_mu<- NA
} else {
  MSE_mu <- get_dist(mu_mat_est, true_mu, dist_type = "mean_sq")
}

#MSE_mu_vec[d]<-get_dist(mu_mat_est, true_mu, dist_type = "mean_sq")

alloc_true_d <- alloc_true[[rid]]
# ARI
if (length(alloc_est) != length(alloc_true_d)) {
  ARI <- NA
} else {
  ARI <- mclust::adjustedRandIndex(alloc_est, alloc_true_d)
}



saveRDS(
  list(
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

# [crodriguezcabrera@boslogin07 ~]$ sbatch Mbeta_Project/my_array_job.sh
# Submitted batch job 142852

