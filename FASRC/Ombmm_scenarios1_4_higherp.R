#=================================================================
# Use simulated data from Cases 1-4 (see excel file with values and see  R/simulate_data_casesv2.R)
#For each case,  we generate 100 datasets  per case to calculate ARI, MSE, Bias, and PIPs
# Programmer: Carmen Rodriguez C.
# Date updated: 4/3/2026
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
#------------Load unlabeled datasets for cases 1-4------------------
datlist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_dpgs1_4_higherp.rds")
alloclist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_dpgs1_4_higherp.rds")
#contains tibble with true parameters


#list of lists of datasets (each scenario has 100 datasets of n=2000)

#list of true allocations per scenario
true_alloc<- list(
  case1 = alloclist[[2]]$case1,
  case2 = alloclist[[2]]$case2,
  case3 = alloclist[[2]]$case3,
  case4 = alloclist[[2]]$case4)

scenariosdat<- alloclist[[1]] #tibble of true parameters (mu, lambda, pi)

# Define all scenario × replicate pairs
scenarios <- names(datlist)  # e.g., "case1": "case4"
jobs <- expand.grid(scen = scenarios, rep = 1:100, stringsAsFactors = FALSE)
cid <- jobs$scen[job_id]
rid <- jobs$rep[job_id]



start_time <- Sys.time()
#changing prior for lambda because we have different values
res<-run_1sim_OmbmmSF(cid, rid, datlist, true_alloc, sceneriosdat,
                      Kmax = 20, niter = 10000, nburnin = 3000,
                      g =  5, h = 10, a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
                      a_rho = 1, b_rho = 20, e0=0.01, priorOnAlpha = "no")


#list of size 3
end_time <- Sys.time()
runtime_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))

res$runtime_sec <- runtime_sec


#Save result 
out_file <- sprintf("/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/Ombmmresult_hp_%s_rep%d.rds", cid, rid)
saveRDS(res, out_file)



# start_time <- Sys.time()
# testing1<-run_1sim_OmbmmSF("case1", 1, datlist, true_alloc, sceneriosdat,
#                  Kmax = 20, niter = 10000, nburnin = 3000,
#                  g =  5, h = 10, a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
#                  a_rho = 1, b_rho = 20, e0=0.1, priorOnAlpha = "no")
# end_time <- Sys.time()
# #
# # # 
# 
# start_time <- Sys.time()
# testing2<-run_1sim_OmbmmSF("case2", 1, datlist, true_alloc, sceneriosdat,
#                            Kmax = 10, niter = 10000, nburnin = 3000,
#                            g =  5, h = 10, a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
#                            a_rho = 1, b_rho = 20, e0=0.1, priorOnAlpha = "no")
# end_time <- Sys.time()
# # 
