#Last updated: 1/19/26
setwd("/n/home03/crodriguezcabrera/bayesmbmm/")

library(tidyverse)
library(mclust)
library(dplyr)
#library(clustvarsel)


datlist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/MockACSUpdated/datlistAll_mockacsdpgs1_5.rds") 
#alloclist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/MockACSUpdated/alloclistAll_mockacsdpgs1_5.rds")

#scenariosdat<- alloclist[[1]] #tibble of true parameters (mu, lambda, pi)
#true_alloc<-alloclist[[2]]
#-------------------------------------------------------------
#-------------------------------------------------------------
#-------------------------------------------------------------
acsmock3<-datlist$MockACS3
#acsmock3_alloc<-true_alloc$MockACS3
#trueparams3<-scenariosdat %>% filter(dgp_id == "MockACS3")

acsmock5<-datlist$MockACS5
#acsmock5_alloc<-true_alloc$MockACS5
#trueparams5<-scenariosdat %>% filter(dgp_id == "MockACS5")

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
mclust_fits3 <- vector("list", nrep)
mclust_fits5 <- vector("list", nrep)

for (r in 1:nrep){
  message("Replication ", r)
  
  mclust_fits3[[r]] <- run_mclust_rep(
    dat = acsmock3[[r]],
    G_grid = 1:20
  )
}

for (r in 1:nrep){
  message("Replication ", r)
  
  mclust_fits5[[r]] <- run_mclust_rep(
    dat = acsmock5[[r]],
    G_grid = 1:20
  )
}

saveRDS(list(mockacs3=mclust_fits3, mockacs5= mclust_fits5),
        file = "/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/mclust_fits_mockacs_new.rds")


