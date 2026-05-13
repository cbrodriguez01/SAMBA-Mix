
library(BayesBinMix)

args <- commandArgs(trailingOnly = TRUE)
rid <- as.integer(args[1])
set.seed(2025 + rid)


#------------Load binarized data: ~/bayesmbmm/bayesbinmix_dataprep.R
case6_binary <- readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/case6_binary.rds")


#---BayesBinMix---
#Using code from by ecbayesbinmix paper
#Function to generate heating vectors-- using method from Altekar 2004
#' @param deltatempvec vector containing multiple deltaT for the heats vector
#' #in the other models smaller delta (not too small) worked best-- we used 0.025
generateHeats<-function(deltatempvec, npchains){
  heatslist<-list()
  for (j in 1:length(deltatempvec)){
    heats<-c()
    for (i in 1:npchains){
      heats[i]<- 1/(1 +  (i-1) * deltatempvec[j])
    }
    #print(heats)
    #Save heats vector in list
    heatslist[[j]]<-heats
  }
  return(heatslist)
  
}

#Model parameters
wd<-"/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/model_comp/"
nChains<- 4
Kmax<-20
m<-5000
burnin<-2000
ClusterPrior <- "poisson"
gamma<- rep(1, Kmax)
deltatemps<-c(0.01,0.025)
heatslist<-generateHeats(deltatemps, nChains)
heatsid<-paste0("deltatemp_", formatC(deltatemps, format = "f", digits = 3))
names(heatslist) <-heatsid


# Run model for this rid
dataset <- case6_binary[[rid]]
for (j in seq_along(heatslist)) {
      heatvec <- heatslist[[j]]
      heatsid <- names(heatslist)[j]
      
      outPrefix <- paste0(
        wd, "case6_rep", rid, "_", heatsid
      )
      
   
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

