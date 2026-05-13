#Examine the sensitivity of Feature Saliency as I change the parameters for rho ~ beta(a,b)
#Remember rho is the probability of inclusion for all variables-- global

args <- commandArgs(trailingOnly = TRUE)
rid <- as.integer(args[1])

set.seed(1000 + rid)

setwd("/n/home03/crodriguezcabrera/bayesmbmm/")
library(nimble)
library(coda)
library(tidyverse)
library(dplyr)
source("R/bayesmbmm_OFMM_SF_v2.R")
source("R/helpers.R")

#updated 4/4/26
datlist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_dpgs1_4_higherp.rds")
  alloclist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_dpgs1_4_higherp.rds")
scenariosdat<- alloclist[[1]]


#DPGs for  C1 AND C2
sim_data<- list(
  case3 = datlist$case3,
  case4 = datlist$case4)

true_params3<- get_true_params(scenariosdat, "case3")
#true_params4<- get_true_params(scenariosdat, "case4")
true.pi<- true_params3$true_pi
#true_mu<- true_params4$true_mu
p<-true_params3$numfeat


#rho ~ beta(a,b) => E(rho) = a/ a+ b: This directly determines how many variables the model expects
#to be salient before seeing data.

# If a=b, E(rho) = 0.5
#if a << b then E(rho) is very small=> induce sparsity
#if a>>b, then E(rho) would be larger and so more variables will be included

# a<-c(0.5,1,2,5,10)
# b<- c(0.5,1,5,10,20)
# hp<-expand_grid(a,b)

#Think about E(rho) instead -- this reduces to 4
# so If E(rho) = 0.05 we'll have very sparse 

rho_mean <- c(0.05, 0.10, 0.25, 0.50)
conc <- 10 #a+b

hp <- tibble(
  a = rho_mean * conc,
  b = (1 - rho_mean) * conc)


#Test this with case 3 first-- C1
case3<-sim_data$case3
alloc_true<- alloclist[[2]][["case3"]]

n_hp <- nrow(hp)

res<- matrix(NA_real_, nrow = n_hp, ncol = 2+p)
colnames(res) <- c("K", "ARI", paste0("PIP", 1:p))

for (i in 1:n_hp){
    a_rho <-as.numeric(hp[i,1])
    b_rho <-as.numeric(hp[i,2])
    message(sprintf("Dataset %f | a_rho=%f | b_rho=%f",
                    rid, a_rho, b_rho))
    
  success<- tryCatch({
  fit <- Ombmm_FeatSaliency(
    x = case3[[rid]],
    K = 20,
    niter = 2000,
    nburnin = 500,
    nchain = 1,
    init_list = NULL, 
    g =  5, h = 3, a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
    a_rho = a_rho, b_rho = b_rho, e0 = 0.01, priorOnAlpha = "no",
    p.true = true.pi, plot_trace = FALSE)
  
    Kmode<-fit$K
    alloc_est<-fit$modeloutput_perchain$clusterMembershipPerMethod[1,] 
    ARI<- mclust::adjustedRandIndex(alloc_est, alloc_true[[rid]])
    PIP = as.vector(fit$FeatureSal_out[[2]])
  
  res[i,]<- c(Kmode,ARI, PIP)
  
  TRUE
  
   },error = function(e){
     message(sprintf(
       "skipping dataset %f | a_rho=%f | b_rho=%f because: %s",
       rid, a_rho, b_rho, conditionMessage(e)
     ))
     
  FALSE
  
   })
   if (success){
     message("Completed successfully.")
   }
  }


out_file <- sprintf(
  "/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/sensitivityrho_hp_case3_rep%d.rds",
  rid
)

saveRDS(
  list(
    rid = rid,
    hyperparams = hp,
    results = res
  ),
  out_file
)

#library(GGally)
#ggpairs(as.data.frame(case4[[100]]))
#plot(as.data.frame(case4[[100]]))

