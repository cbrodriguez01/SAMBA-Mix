library(BayesBinMix)

args <- commandArgs(trailingOnly = TRUE)
rid <- as.integer(args[1])
set.seed(2025 + rid)


#------------Load binarized data: ~/bayesmbmm/bayesbinmix_dataprep.R
sim_data_binary<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/sim_data_binary_mockacs3_5.rds")


#---BayesBinMix---

#Model parameters
wd<-"/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/model_comp/"
nChains<- 4
Kmax<-20
m<-5000
burnin<-2000
ClusterPrior <- "poisson"
gamma<- rep(1, Kmax)
#deltatemp<-0.01# using approach from ecbayesbinmix paper: generateHeats(0.01, nChains) 
heatvec<-c(1, 0.9900990,0.9803922, 0.9708738)

# Run model for this rid

for (j in names(sim_data_binary)) {
  
  dataset<-sim_data_binary[[j]][[rid]]
  
  outPrefix <- paste0(
    wd, "BayesBinMix_ACSMOCK_",j, "_rep", rid)
  
  
  res <- coupledMetropolis(
    Kmax = Kmax,
    nChains = nChains,
    heats = heatvec,
    binaryData = dataset,
    outPrefix = outPrefix,
    ClusterPrior = ClusterPrior,
    m = m,
    gamma = gamma,
    burn = burnin
  )
  
  # Save raw output only
  saveRDS(
    res,
    file = paste0(outPrefix, ".rds")
  )
}
