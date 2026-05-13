#=================================================================
# Use simulated data from ACSmock (see ~/bayesmbmm/data/simulate_acs_mockdatav2.R)
#For each case,  we generate 100 datasets 
# Programmer: Carmen Rodriguez C.
# Date updated: 3.10.26
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


datlist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/MockACSUpdated/datlistAll_mockacsdpgs1_5.rds") 
alloclist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/MockACSUpdated/alloclistAll_mockacsdpgs1_5.rds")

scenariosdat<- alloclist[[1]] #tibble of true parameters (mu, lambda, pi)

# Define all scenario × replicate pairs
scenarios <- names(datlist)  
jobs <- expand.grid(scen = scenarios, rep = 1:100, stringsAsFactors = FALSE)
cid <- jobs$scen[job_id]
rid <- jobs$rep[job_id]
casedat<-datlist[[cid]]
alloc_true <- alloclist[[2]][[cid]]
true_params<- get_true_params(scenariosdat, cid)
p<-true_params$numfeat


#Save result 
out_file <- sprintf("/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/SAMBOACSmockFinal_%s_rep%d.rds", cid, rid)

start_time <- Sys.time()
# Fit model
fit <- tryCatch(
  {
    Ombmm_FeatSaliency(
      x = casedat[[rid]],
      K = 20,
      niter = 20000, 
      nburnin = 5000,
      nchain = 1,
      init_list = NULL, 
      g =  5, h = 3, 
      a_mu0 = 1, b_mu0 = 1, 
      a_mu1 = 1, b_mu1 = 1,
      a_rho = 1, b_rho = 1, 
      e0=1, priorOnAlpha = "no",
      p.true = true_params$true_pi,
      plot_trace = FALSE
    )
    
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

alloc_true_d <- alloc_true[[rid]]
# ARI
if (length(alloc_est) != length(alloc_true_d)) {
  ARI <- NA
} else {
  ARI <- mclust::adjustedRandIndex(alloc_est, alloc_true_d)
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


#Submitted batch job 66568210

