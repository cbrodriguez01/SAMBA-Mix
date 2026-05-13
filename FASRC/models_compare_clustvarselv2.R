args <- commandArgs(trailingOnly = TRUE)
rid <- as.integer(args[1])   # replication index (1–100)


cat("Running replication", rid, "\n")


#Last updated: 3/18/26
setwd("/n/home03/crodriguezcabrera/bayesmbmm/")
#https://cran.r-project.org/web/packages/clustvarsel/vignettes/clustvarsel.html

library(tidyverse)
library(mclust)
library(dplyr)
library(clustvarsel)

set.seed(2025 + rid)

datlist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/MockACSUpdated/datlistAll_mockacsdpgs1_5.rds") 

acsmock3<-datlist$MockACS3[[rid]]
acsmock5<-datlist$MockACS5[[rid]]


#-------------------------------------------------------------
#-------------------------------------------------------------
#-------------------------------------------------------------

run_clustvarsel_rep<-function(dat,G_grid= 2:20, drt){
  
  x<- apply(dat, 2, qlogis)
  
  # Fit model
  t0 <- proc.time()
  mc <- clustvarsel(x, G = G_grid, direction = drt)
  
  runtime <- (proc.time() - t0)["elapsed"]
  
  
  modelout<- mc$model #fitted model with the selected variables
  
  list(
    subset.        =  mc$subset,  #variables selected
    classification = modelout$classification,
    G_hat          = modelout$G,
    modelName      = modelout$modelName,
    BIC            = modelout$BIC,
    parameters     = modelout$parameters,
    runtime        = runtime
  )
}


fit3 <- run_clustvarsel_rep(
  dat    = acsmock3,
  G_grid = 1:20,
  drt    = "backward"
)

fit5 <- run_clustvarsel_rep(
  dat    = acsmock5,
  G_grid = 1:20,
  drt    = "backward"
)

outdir <- "/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs"
outfile <- file.path(outdir,
                     sprintf("clustvarsel_mockacs_rep%d.rds", rid)
)

saveRDS(
  list(mockacs3=fit3, mockacs5= fit5),
  outfile
)

cat("Saved:", outfile, "\n")

