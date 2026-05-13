#make binary version of case 6 data--only load once
library(tidyverse)
library(dplyr)


#UPDATED: 3/18/26
#USED OLD DATASETS IN SIMS 
# Load new datasets
datlist <- readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/MockACSUpdated/datlistAll_mockacsdpgs1_5.rds") 
alloclist <- readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/MockACSUpdated/alloclistAll_mockacsdpgs1_5.rds")

# True parameter tibble
scenariosdat <- alloclist[[1]]

# Extract only the datasets you want
acsmock3 <- datlist$MockACS3
acsmock5 <- datlist$MockACS5

# Build sim_data in the same structure as before
sim_data <- list(
  MockACS3 = acsmock3,
  MockACS5 = acsmock5
)

#For BAYESBINMIX -- we do binary (median)

binarize_case <- function(case_list) {
  medians_list <- lapply(case_list, function(x) {
    apply(x, 2, median)
  })
  
  bindat <- lapply(seq_along(case_list), function(i) {
    x   <- case_list[[i]]
    med <- medians_list[[i]]
    sweep(x, 2, med, FUN = ">=") * 1
  })
  
  return(bindat)
}

sim_data_binary <- lapply(sim_data, binarize_case)


saveRDS(
  sim_data_binary,
  file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/sim_data_binary_mockacs3_5.rds"
)




#############################################################################################
#ADDED 4/3/26
#############################################################################################
#For VICATMIX -- we do tertiles
datlist7<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_case7.rds") 
datlist5<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_case5.rds")
sim_data<- list(case7 = datlist7, case5= datlist5)


compute_global_tertiles <- function(case_list) {
  x_all <- do.call(rbind, case_list)
  
  apply(x_all, 2, quantile, probs = c(1/3, 2/3), na.rm = TRUE)
}


apply_tertiles <- function(case_list, cuts) {
  
  lapply(case_list, function(x) {
    
    out <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
    
    for (j in 1:ncol(x)) {
      out[, j] <- cut(
        x[, j],
        breaks = c(-Inf, cuts[1, j], cuts[2, j], Inf),
        labels = c(0, 1, 2),
        include.lowest = TRUE
      )
    }
    
    # VICatMix expects integers
    out <- apply(out, 2, function(col) as.integer(as.character(col)))
    
    return(out)
  })
}

sim_data_tertile <- lapply(sim_data, function(case_list) {
  cuts <- compute_global_tertiles(case_list)
  apply_tertiles(case_list, cuts)
})



saveRDS(
  sim_data_tertile,
  file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/sim_data_tertile_case5_case7.rds"
)





















#------------Load unlabeled datasets for cases 1-5------------------
#datlist6<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_case6.rds")
datlist7<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_case7.rds") 
datlist5<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_case5.rds")

alloclist6<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case6.rds")
alloclist7<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case7.rds")
alloclist5<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case5.rds")

#case6 = datlist6[[1]]
#case6_alloc<-alloclist6[[2]]

scenariosdat<- alloclist6[[1]] #tibble of true parameters (mu, lambda, pi)
trueparams7<-scenariosdat %>% filter(cid == "case7")
trueparams5<-scenariosdat %>% filter(cid == "case5")
sim_data<- list(case7 = datlist7, case5= datlist5)
true_alloc<- list(case7=alloclist7, case5=alloclist5)

#---BayesBinMix---
#https://www.reddit.com/r/statistics/comments/1gs225j/q_how_to_define_a_cut_off_value_to_create_binary/
# Transform all  all variables to binary
# var_medians_list<- lapply(case6,function(i) apply(i, 2, median))
# 
# #Create the binary variables, 1 if above or equal to median and 0 if below 
# case6_binary <- lapply(seq_along(case6), function(i) {
#   x <- case6[[i]]
#   med <- var_medians_list[[i]]
#   
#   sweep(x, 2, med, FUN = ">=") * 1
# })


binarize_case <- function(case_list) {
medians_list <- lapply(case_list, function(x) {
    apply(x, 2, median)
  })
  
bindat <- lapply(seq_along(case_list), function(i) {
    x   <- case_list[[i]]
    med <- medians_list[[i]]
    sweep(x, 2, med, FUN = ">=") * 1
  })
  
  return(bindat)
}

sim_data_binary <- lapply(sim_data, binarize_case)


saveRDS(
  sim_data_binary,
  file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/sim_data_binary_case5_case7.rds"
)

#saveRDS(case6_binary, file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/case6_binary.rds")
