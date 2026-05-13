library(tidyverse)
library(nimble)
source("~/bayesmbmm/R/helpers.R")
source("~/bayesmbmm/R/bayesmbmm_OFMM_SF_v2.R")

args <- commandArgs(trailingOnly = TRUE)
job_id <- as.numeric(args[1])

set.seed(2026 + job_id)



acs<-readRDS(file = "~/AHRQData/acs19.rds")
acs <- as.matrix(acs)
storage.mode(acs) <- "double"

configs <- list(
  list(name= "conf1", e0=0.001, alpha_prior="no", a_rho=0.5,   b_rho=15),
  list(name= "conf2", e0=0.01,  alpha_prior="no", a_rho=1,   b_rho=10),
  list(name= "conf3", e0=0.1,   alpha_prior="gam_1_2", a_rho=1, b_rho=10),
  list(name= "conf4", e0=0.001, alpha_prior="no", a_rho=1,   b_rho=1)
)

current_config <- configs[[job_id]]
message(paste("Running config:", current_config$name))



start_time <- Sys.time()
fit <- Ombmm_FeatSaliency(
  x = acs,
  K = 30,
  niter = 15000,
  nburnin = 5000,
  nchain = 1,
  init_list = NULL,
  g = 10, h = 2,      # Tighter clusters
  a_mu0 = 2, b_mu0 = 2, 
  a_mu1 = 2, b_mu1 = 2, 
  a_rho = current_config$a_rho, 
  b_rho = current_config$b_rho,
  e0 = current_config$e0,
  priorOnAlpha = current_config$alpha_prior,
  plot_trace = FALSE
)

end_time <- Sys.time()
runtime_sec <- as.numeric(difftime(end_time, start_time, units = "secs")) 

# --- SAVE RESULTS ---
base_dir <- "/n/netscratch/stephenson_lab/Lab/crodriguez/AHRQ_application/"
out_file <- sprintf("%sacs19_MA_sambo_%s.rds", base_dir, current_config$name)

saveRDS(list(
  model_fit = fit, 
  runtime = runtime_sec,
  config = current_config
), out_file)



#We need priors that
#1. respect skweness
#2. allow cluster variation
#3. dont over-regularize lambda
#4. encourage some sparsity in features but not too much


#Priors for mu
#these are mostly small proportions so we need to make the prior can explore these spaces
#lambda
# acs vars have different dispersion levels

