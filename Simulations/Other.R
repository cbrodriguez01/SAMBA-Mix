
#To generate plots for the supplementary materials
# Load libraries

setwd("/n/home03/crodriguezcabrera/bayesmbmm/")
library(nimble)
library(coda)
library(tictoc)
library(tidyverse)
library(dplyr)
source("R/helpers.R")
source("R/bayesmbmm_OFMM_SF_v2.R")

#------------Load unlabeled datasets for cases 1-4------------------
datlist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_dpgs1_4.rds")
alloclist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_dpgs1_4.rds")
#contains tibble with true parameters

#Randomly select one dataset from each case
set.seed(2026)
rid<- sample(100, size = 1)
cid<-names(datlist)

sim_data<- lapply(datlist, function(x) x[[rid]])


#Save fit
res <- vector("list", length(sim_data))
names(res) <- names(sim_data)

for (i in seq_along(sim_data)){
  fit<-Ombmm_FeatSaliency(
    x = sim_data[[i]],
    K = 10,
    niter = 10000, 
    nburnin = 3000,
    nchain = 1,
    init_list = NULL, 
    g =  5, h = 3, a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
    a_rho = 1, b_rho = 20, e0 = 0.01,
    priorOnAlpha = "no",
    p.true = NULL, plot_trace = TRUE)
  
  res[[i]]<-fit
  
}
  
saveRDS(res, file= "~/Mbeta_Project/fortraceplots.rds")
  



