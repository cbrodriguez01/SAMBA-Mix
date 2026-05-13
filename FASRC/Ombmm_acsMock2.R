#=================================================================
#Goal: Having underfitting issues... and so we want to test
# Is underfitting driven by e0?
#   Is it amplified by variable selection?
#   Do they interact?
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

datlist5<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_case5.rds")
sim_dat<-datlist5[[24]]

alloclist6<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case6.rds")
scenariosdat<- alloclist6[[1]] #tibble of true parameters (mu, lambda, pi)
scenariosdat<-scenariosdat %>% rename(dgp_id = cid)
true_params<- get_true_params(scenariosdat, "case5")
true.pi<- true_params$true_pi


# Define all scenario × replicate pairs
jobs <- expand.grid(
  e0 = c(0.01, 0.05, 0.1),
  rho_a = c(1, 5, 10),
  rho_b = 1
)

this_hp <- jobs[job_id, ]
e0 <- this_hp$e0
a_rho <-this_hp$rho_a
b_rho <-this_hp$rho_b

start_time <- Sys.time()
fit <- Ombmm_FeatSaliency(
  x = sim_dat,
  K = 20,
  niter = 12000,
  nburnin = 4000,
  nchain = 1,
  init_list = NULL, 
  g =  5, h = 3, a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
  a_rho = a_rho, b_rho = b_rho, e0 = e0, priorOnAlpha = "no",
  p.true = true.pi, plot_trace = FALSE)
end_time <- Sys.time()
runtime_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))

out <- list(
  fit = fit,
  e0 = e0,
  a_rho = a_rho,
  b_rho = b_rho,
  runtime_sec = runtime_sec
)

#Save result 
out_file <- sprintf(
  "/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/fit_e0_%s_rho_%s_%s_job%d.rds",
  e0, a_rho, b_rho, job_id
)
saveRDS(out, out_file)


#------First we look at model parameters for all
files_path<-"/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs"

test<-readRDS("/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/fit_e0_0.01_rho_1_1_job1.rds")
out<-test$fit



head(out$nk_iter)
colMeans(out$nk_iter)





fit <- Ombmm_FeatSaliency(
  x = sim_dat,
  K = 5,
  niter = 12000,
  nburnin = 4000,
  nchain = 1,
  init_list = NULL, 
  g =  5, h = 3, a_mu0 = 1, b_mu0 = 1,a_mu1 = 1, b_mu1 = 1,
  a_rho = a_rho, b_rho = b_rho, e0 = e0, priorOnAlpha = "no",
  p.true = true.pi, plot_trace = FALSE)




scenariosdat[3,]$lambda_vec

scenariosdat[3,]$mu
