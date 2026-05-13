library(tidyverse)
library(nimble)
source("~/bayesmbmm/R/helpers.R")
source("~/bayesmbmm/R/bayesmbmm_OFMM_SF_v2.R")

# 
# args <- commandArgs(trailingOnly = TRUE)
# job_id <- as.numeric(args[1])

set.seed(2030)


acs<-readRDS(file = "~/AHRQData/acs19.rds")
acs <- as.matrix(acs)
storage.mode(acs) <- "double"

start_time <- Sys.time()
fit <- Ombmm_FeatSaliency(
  x = acs,
  K = 30,
  niter = 20000,
  nburnin = 10000,
  nchain = 1,
  init_list = NULL,
  g = 2, h = 1,     # Tighter clusters
  a_mu0 = 2, b_mu0 = 2, 
  a_mu1 = 2, b_mu1 = 2, 
  a_rho = 2, 
  b_rho = 5,
  e0 = 0.001,
  priorOnAlpha = "no",
  plot_trace = FALSE
)

end_time <- Sys.time()
runtime_sec <- as.numeric(difftime(end_time, start_time, units = "secs")) 


out_file <- sprintf("/n/netscratch/stephenson_lab/Lab/crodriguez/AHRQ_application/acs19_sambo_4.8.26.rds")

saveRDS(list(fit, runtime_sec), out_file)



res <- readRDS("/n/netscratch/stephenson_lab/Lab/crodriguez/AHRQ_application/acs19_sambo_4.8.26.rds")
fit <- res[[1]]
runtime <- res[[2]]/3600 # 10 hours
Kmode<-fit$K
fit$FeatureSal_out[[2]]
# s[1]  s[2]  s[3]  s[4]  s[5]  s[6]  s[7]  s[8]  s[9] s[10] s[11] s[12] s[13] s[14] 
# 0     0     0     1     0     1     0     1     0     1     1     0     0     1 

names(acs)
# 
# > names(acs)
# [1] "Femalehousehold_2019_P"    "less_than_hs_p"            "Crowding_housing_2019"    
# [4] "working_class_2019"        "UnemployementP_2019"       "RenterOccupiedUnitP_2019" 
# [7] "No_vehicle_2019"           "Lackplumbing_2019"         "NonHispanicBlack_2019"    
# [10] "NonHispanicAsian_2019"     "Hispanic_or_Latino_2019"   "lang_home_EN_notwell_2019"
# [13] "SNAP_2019_P"               "income_scaled"            
