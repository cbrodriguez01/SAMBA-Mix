#Last updated: 1/19/26
setwd("/n/home03/crodriguezcabrera/bayesmbmm/")

library(tidyverse)
library(mclust)
library(dplyr)
#library(clustvarsel)

set.seed(2025)
#-Note: case 7 contains 20 varaibles where 12 are important and 8 are noise; case5 contains only the 12 important variables
#datlist6<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_case6.rds")
alloclist6<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case6.rds")

datlist7<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_case7.rds") 
datlist5<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_case5.rds")


alloclist7<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case7.rds")
alloclist5<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case5.rds")

sim_data<- list(case7 = datlist7, case5= datlist5)
#list of true allocations per scenario
true_alloc<- list(case7=alloclist7, case5=alloclist5)
scenariosdat<- alloclist6[[1]] #tibble of true parameters (mu, lambda, pi)


#-------------------------------------------------------------
#-------------------------------------------------------------
#-------------------------------------------------------------
case5<-sim_data$case5
case5_alloc<-true_alloc$case5
trueparams5<-scenariosdat %>% filter(cid == "case5")

case7<-sim_data$case7
case7_alloc<-true_alloc$case7
trueparams7<-scenariosdat %>% filter(cid == "case7")

#-----Gaussian Mixture Model-----
#Data needs to be transformed before using mclust; we can also try clustvarsel
# 
# logit <- function(x, eps = 1e-6) {
#   x <- pmin(pmax(x, eps), 1 - eps)
#   log(x / (1 - x))
# }
# 
# inv_logit <- function(z) exp(z) / (1 + exp(z))

#Note: We are using EEE, and VVV for the covariance structure--see book
#EEE assumes identical cluster structure and VVV allows for varying 
#https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html#clustering
#Book: https://mclust-org.github.io/mclust-book/index.html

run_mclust_rep<-function(dat,G_grid= 2:20){
  
  x<- apply(dat, 2, qlogis)
  
  # Fit mclust
  t0 <- proc.time()
  mc <- Mclust(
    x,
    G = G_grid,
    modelNames = c("EEE", "VVV")
  )
  runtime <- (proc.time() - t0)["elapsed"]
  
  
  list(
    classification = mc$classification,
    G_hat          = mc$G,
    modelName      = mc$modelName,
    BIC            = mc$BIC,
    parameters     = mc$parameters,
    runtime        = runtime
  )
}

#i could do this better but for now this is okay!
nrep <- 100
mclust_fits5 <- vector("list", nrep)
mclust_fits7 <- vector("list", nrep)

for (r in 1:nrep){
  message("Replication ", r)
  
  mclust_fits5[[r]] <- run_mclust_rep(
    dat = case5[[r]],
    G_grid = 1:20
  )
}

for (r in 1:nrep){
  message("Replication ", r)
  
  mclust_fits7[[r]] <- run_mclust_rep(
    dat = case7[[r]],
    G_grid = 1:20
  )
}

saveRDS(list(case5=mclust_fits5, case7= mclust_fits7),
  file = "/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/mclust_fits_cases5and7.rds")


