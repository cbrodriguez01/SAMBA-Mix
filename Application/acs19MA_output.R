setwd("~/bayesmbmm/results")
library(nimble)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(coda)
library(bayesplot)
library(flextable)
library(officer)
library(writexl)
library(tigris) # To get MA shapefiles
library(sf)

library(tidyverse)
library(knitr)
library(kableExtra)

base_dir <- "/n/netscratch/stephenson_lab/Lab/crodriguez/AHRQ_application/"
configs <- c("conf1", "conf2", "conf3", "conf4")

# load and extract 
summary_list <- list()

for(cfg in configs){
  file_path <- paste0(base_dir, "acs19_MA_sambo_", cfg, ".rds")
  
  if(file.exists(file_path)){
    res <- readRDS(file_path)
    fit <- res$model_fit
    

    Kmode <- fit$K
    
    # Feature Saliency (Count variables with PIP > 0.5)
    pip<-as.numeric(fit$FeatureSal_out[[2]])
    pip_count<-sum(pip > 0.5)
    
    summary_list[[cfg]] <- data.frame(
      Configuration = cfg,
      e0 = res$config$e0,
      Kmode = Kmode,
      Salient_Features = pip_count,
      Runtime_Min = round(res$runtime / 60, 2)
    )
  }
}

final_comparison <- do.call(rbind, summary_list)


file<-paste0(base_dir, "acs19_MA_sambo_", "conf1", ".rds")
res <- readRDS(file)
fit <- res$model_fit
fit$FeatureSal_out[[2]]


#conf3
# [1] "/n/netscratch/stephenson_lab/Lab/crodriguez/AHRQ_application/acs19_MA_sambo_conf3.rds"
# > res <- readRDS(file)
# > fit <- res$model_fit
# > fit$FeatureSal_out[[2]]
# s[1]  s[2]  s[3]  s[4]  s[5]  s[6]  s[7]  s[8]  s[9] s[10] s[11] s[12] s[13] s[14] 
# 0     1     0     1     0     1     1     0     0     0     1     1     1     1 
# 
# > names(acs)
# [1] "Femalehousehold_2019_P"    "less_than_hs_p"            "Crowding_housing_2019"     "working_class_2019"       
# [5] "UnemployementP_2019"       "RenterOccupiedUnitP_2019"  "No_vehicle_2019"           "Lackplumbing_2019"        
# [9] "NonHispanicBlack_2019"     "NonHispanicAsian_2019"     "Hispanic_or_Latino_2019"   "lang_home_EN_notwell_2019"
# [13] "SNAP_2019_P"               "income_scaled"            
