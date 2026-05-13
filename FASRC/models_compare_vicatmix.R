# VICatMix variable selection for binarized simulation data
#https://cran.r-project.org/web/packages/VICatMix/VICatMix.pdf
#See paper for more info
library(VICatMix)
args <- commandArgs(trailingOnly = TRUE)
rid <- as.integer(args[1])
set.seed(2025 + rid)
#https://cran.r-project.org/web/packages/VICatMix/VICatMix.pdf

#------------Load binarized data: ~/bayesmbmm/bayesbinmix_dataprep.R
sim_data_tertile<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/sim_data_tertile_case5_case7.rds")
wd<-"/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/"



for (j in names(sim_data_tertile)) {
  dataset<-sim_data_tertile[[j]][[rid]]
  res <- runVICatMixVarSel(
    data = dataset,
    K = 20,
    alpha = 0.01)
  if (is.null(res) || all(is.na(res$labels))) {
    warning(paste("VICatMix failed for case", j, "rep", rid))
    next}
  #   
  # vicatmix_out[[j]][[rid]]<- list(classification = res$labels,
  #             clustersnum = res$Cl, #number of clusters in every iteration
  #             pip = res$model$c #A P-length vector of expected values for the variable selection parameter, gamma.
  # )
  
  saveRDS(
    list(
      classification = res$labels,
      clustersnum = res$Cl,
      pip = res$model$c
    ),
    file = paste0(wd, "VICatMix_cases5_7_tertiles", j, "_rep", rid, ".rds")
  )
}
  



# save one file per case 
# for (j in names(vicatmix_out)) {
#   saveRDS(
#     vicatmix_out[[j]],
#     file = paste0(wd, "VICatMix_", j, "_allreps.rds")
#   )
# }

