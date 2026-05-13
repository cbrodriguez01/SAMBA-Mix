#Updated to add more repetitions (2/9/26)
#Updated 4/4/26 to change datasets
#Look through priors e0 and varying number of Kmax using dataset case 4 == scenario_id C2, 

# Load libraries
setwd("/n/home03/crodriguezcabrera/bayesmbmm/")
library(nimble)
library(coda)
library(tidyverse)
library(dplyr)
source("R/bayesmbmm_OFMM_SF_v2.R")
source("R/helpers.R")

e0_vec <- c(0.1, 0.01, 0.005)
priorOnAlpha_vec<- c("gam_05_05", "gam_1_2")
Kmax_vec<-c(10,30)

# Case 1: no prior on alpha
hp_no_alpha <- expand.grid(
  Kmax = Kmax_vec,
  e0 = e0_vec,
  priorOnAlpha = "no",
  stringsAsFactors = FALSE
)

# Case 2: prior on alpha 
hp_with_alpha <- expand.grid(
  Kmax = Kmax_vec,
  e0 =NA,   
  priorOnAlpha = priorOnAlpha_vec,
  stringsAsFactors = FALSE
)

hp <- dplyr::bind_rows(hp_no_alpha, hp_with_alpha)
#nrow(hp); 10
hp_rid <- expand.grid(
  rep = 1:50,
  seq_len(nrow(hp)),
  stringsAsFactors = FALSE
)
hp_rid <- cbind(
  hp_rid["rep"],
  hp[hp_rid[[2]], ]
)
rownames(hp_rid) <- NULL


# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
hpid<-as.numeric(args[1])

this_hp <- hp_rid[hpid, ]
Kmax <- this_hp$Kmax
e0 <- this_hp$e0
priorOnAlpha <- this_hp$priorOnAlpha
rid <- this_hp$rep


set.seed(1000 + hpid)

datlist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_dpgs1_4_higherp.rds")
alloclist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_dpgs1_4_higherp.rds")
scenariosdat<- alloclist[[1]]

sim_data<- datlist$case4
true_params<- get_true_params(scenariosdat, "case4")
true.pi<- true_params$true_pi
true_mu<- true_params$true_mu
p<-true_params$numfeat
alloc_true<- alloclist[[2]]$case4


#Fit model
    fit <- tryCatch({
      
      Ombmm_FeatSaliency(
      x = sim_data[[rid]],
      K = Kmax,
      niter = 2000,
      nburnin = 500,
      nchain = 1,
      init_list = NULL, 
      g =  5, h = 3, a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
      a_rho = 1, b_rho = 20, e0 = if(priorOnAlpha == "no") e0 else 0.01,
      priorOnAlpha = priorOnAlpha,
      p.true = NULL, plot_trace = FALSE)} ,

      error = function(e){
       message(sprintf("Job %d | Dataset %d failed: %s",hpid, rid, e$message
       ))
       return(NULL)
     })
    
    
    if (is.null(fit)) {
      Kmode_vec[d] <- 1
      MSE_mu_vec[d] <- NA
      ARI_vec[d] <- NA
      next
    }
    
    Kmode<-fit$K
    alloc_est<-fit$modeloutput_perchain$clusterMembershipPerMethod[1,] #Allocation after ECR-ITERATIVE-1 method for label switching
    Kmode_vec<- Kmode
    sSum <- samplesSummary(fit$modeloutput_perchain$parameters.ECR.mcmc)
    mu_mat_est <- sapply(1:p, function(j) {
      as.numeric(sSum[paste0("mu[", 1:Kmode, ", ", j, "]"), 1])})
    
    if (nrow(mu_mat_est) != nrow(true_mu) ||
        ncol(mu_mat_est) != ncol(true_mu)) {
      MSE_mu_vec <- NA
    } else {
      MSE_mu_vec <- get_dist(mu_mat_est, true_mu, dist_type = "mean_sq")
    }
    
     #MSE_mu_vec[d]<-get_dist(mu_mat_est, true_mu, dist_type = "mean_sq")
     
     alloc_true_d <- alloc_true[[rid]]
     # ARI
     if (length(alloc_est) != length(alloc_true_d)) {
       ARI_vec <- NA
     } else {
       ARI_vec <- mclust::adjustedRandIndex(alloc_est, alloc_true_d)
     }
  



out_file <- sprintf(
  "/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/Sensitivitye0_hp_case4_K%d_e0_%s_%s_rid_%s.rds",
  Kmax, 
  ifelse(is.na(e0), "NA", format(e0, scientific = TRUE)),
  priorOnAlpha, rid)


dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

ok <- tryCatch({
  saveRDS(
    list(
      hyperparams = this_hp,
      Kmode = Kmode_vec,
      MSE_mu = MSE_mu_vec,
      ARIest = ARI_vec
    ),
    out_file
  )
  TRUE
}, error = function(e) {
  message(sprintf("HP %d failed during saveRDS: %s", hpid, e$message))
  FALSE
})

if (!ok) quit(status = 1)
quit(status = 0)




