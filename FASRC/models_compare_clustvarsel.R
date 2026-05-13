args <- commandArgs(trailingOnly = TRUE)
rid <- as.integer(args[1])   # replication index (1–100)

         
cat("Running replication", rid, "\n")


#Last updated: 1/19/26
setwd("/n/home03/crodriguezcabrera/bayesmbmm/")
#https://cran.r-project.org/web/packages/clustvarsel/vignettes/clustvarsel.html

library(tidyverse)
library(mclust)
library(dplyr)
library(clustvarsel)

set.seed(2025 + rid)

#-Note: case 7 contains 20 variables where 12 are important and 8 are noise; case5 contains only the 12 important variables
datlist7<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_case7.rds") 
datlist5<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_case5.rds")


case5<-datlist5[[rid]]
case7<- datlist7[[rid]]


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


fit5 <- run_clustvarsel_rep(
  dat    = case5,
  G_grid = 1:20,
  drt    = "backward"
)

fit7 <- run_clustvarsel_rep(
  dat    = case7,
  G_grid = 1:20,
  drt    = "backward"
)

outdir <- "/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs"
outfile <- file.path(outdir,
                     sprintf("clustvarsel_rep%d.rds", rid)
)

saveRDS(
  list(case5 = fit5, case7 = fit7),
  outfile
)

cat("Saved:", outfile, "\n")







# outdir <- "/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/clustvarsel"
# 
# files <- list.files(outdir, pattern = "clustvarsel_rep", full.names = TRUE)
# files <- files[order(files)]
# 
# fits <- lapply(files, readRDS)
# 
# clustvarsel_fits5 <- lapply(fits, `[[`, "case5")
# clustvarsel_fits7 <- lapply(fits, `[[`, "case7")
# 
# saveRDS(
#   list(case5 = clustvarsel_fits5, case7 = clustvarsel_fits7),
#   "/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/clustvarsel_fits_cases5and7.rds"
# )
# 
# 

