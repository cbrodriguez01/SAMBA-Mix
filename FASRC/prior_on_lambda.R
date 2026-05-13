#Prior sensitivity analysis for concentration parameter lambda
#Using data from scen: case 3 lambda = (25,35) and case 4 lambda =(10,5)
job_id <- as.integer(commandArgs(trailingOnly = TRUE)[1])
set.seed(1000 + job_id)

setwd("/n/home03/crodriguezcabrera/bayesmbmm/")
library(nimble)
library(coda)
library(tidyverse)
library(dplyr)
source("R/bayesmbmm_OFMM_SF_v2.R")
source("R/helpers.R")

#Revised on 4/4/26 to use the p=5 and p=10 datasets
datlist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_dpgs1_4_higherp.rds")
alloclist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_dpgs1_4_higherp.rds")
scenariosdat<- alloclist[[1]]

hp <- tibble(
  g = c(5,1,20, 5, 10,1, 5),
  h = c(1,5, 1, 10, 3, 20, 3))
 
jobs<- expand.grid(
   scen = c("case3", "case4"),
   rep = 1:50,
   hp_id = seq_len(nrow(hp)))

jobs <- cbind(
  jobs[c("scen","rep")],
  hp[jobs$hp_id, ]
)
this_hp <- jobs[job_id, ]

cid<- as.character(this_hp$scen)
casedat<-datlist[[cid]]
rid <- this_hp$rep
g<-this_hp$g
h<-this_hp$h

true_params<- get_true_params(scenariosdat, cid)
p<-true_params$numfeat
alloc_true <- alloclist[[2]][[cid]]


out_file <- sprintf(
  "/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/Sensitivitylambda_hp_%s_jobid_%d_rep_%02d.rds",
  cid, 
  job_id,
  rid)

# Fit model
fit <- tryCatch(
  {
    Ombmm_FeatSaliency(
      x = casedat[[rid]],
      K = 10,
      niter = 5000,
      nburnin = 2000,
      nchain = 1,
      init_list = NULL, 
      g = g, 
      h = h, 
      a_mu0 = 1, b_mu0 = 1,
      a_mu1 = 2, b_mu1 = 2,
      a_rho = 1, b_rho = 20,
      e0 = 0.01,
      priorOnAlpha = "no",
      p.true = true_params$true_pi,
      plot_trace = FALSE
    )
  },
  error = function(e) {
    message(sprintf("Job %d | Dataset %d failed: %s", job_id, rid, e$message))
    return(NULL)
  }
)

if (is.null(fit)) {
  saveRDS(
    list(
      hyperparams = this_hp,
      Kmode = NA,
      mse_lambda = NA,
      mab_lambda = NA,
      MSE_mu = NA,
      ARI = NA,
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
        mab_lambda <- NA
      } else {
        mse_lambda<-get_dist(lambda_est, true_params$true_lambda, dist_type = "mean_sq")
        mab_lambda<-get_dist(lambda_est, true_params$true_lambda, dist_type = "mean_abs")
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
          hyperparams = this_hp,
          job_id = job_id,
          Kmode = Kmode,
          lambda_est = lambda_est, 
          mse_lambda = mse_lambda, 
          mab_lambda = mab_lambda,
          MSE_mu = MSE_mu,
          ARI = ARI),
        out_file
      )

 



# 
# fit <- Ombmm_FeatSaliency(
#   x = sim_data$case3[[1]],
#   K = 10,
#   niter = 2000,
#   nburnin = 500,
#   nchain = 1,
#   init_list = NULL, 
#   g = 1, h = 40 , a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
#   a_rho = 1, b_rho = 20, e0 =  0.01,
#   priorOnAlpha = "no",
#   p.true = true_params3$true_pi, plot_trace = TRUE)
# 
# 
# fit1 <- Ombmm_FeatSaliency(
#   x = sim_data$case4[[1]],
#   K = 10,
#   niter = 2000,
#   nburnin = 500,
#   nchain = 1,
#   init_list = NULL, 
#   g = 1, h = 40 , a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
#   a_rho = 1, b_rho = 20, e0 =  0.01,
#   priorOnAlpha = "no",
#   p.true = true_params4$true_pi, plot_trace = TRUE)


# ## plot trace of lambda
# library(bayesplot) # for color_scheme_set
# samplesSummary(fit$modeloutput_perchain$parameters.ECR.mcmc)
# samplesSummary(fit1$modeloutput_perchain$parameters.ECR.mcmc)
# 
# 
# s<-as.matrix(fit1$modeloutput_perchain$parameters.ECR.mcmc)
# 
# 
# lambda_samples<-s[,paste0("lambda[", 1:2, "]")]
# color_scheme_set("brightblue")
# mcmc_trace(lambda_samples, 
#            pars = paste0("lambda[", 1:2, "]"),
#            facet_args = list(ncol = 1)) +
#   ggtitle("Trace plots for lambda")

