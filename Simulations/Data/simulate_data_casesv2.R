#===================================================================
# Simulate data cases updated
# Programmer: Carmen Rodriguez C.
# Date created: 2025/12/29
#Updated: 4/03/29 -->line 116
#===================================================================
#setwd("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/")

library(tidyverse)
library(ggplot2)
library(mclust)
library(psych)
library(GGally)
library(MASS)
library(corrplot)
source("R/helpers.R")

set.seed(2025)
#-----------DATA CASES-----------------
mu_p3 <-matrix(c(0.2, 0.5, 0.8,   # component 1
                 0.7, 0.3, 0.4),  # component 2
               nrow = 2, byrow = TRUE)


#All variables informative
mu_p5_all<-matrix(c(0.2, 0.5, 0.8,0.4,0.5, 
                    0.7, 0.3, 0.4, 0.5, 0.9),
                  nrow = 2, byrow = TRUE)

#Only 3 variables informative
mu_p5_three<-matrix(c(0.2,0.5,0.8,0.4,  0.6, 
                      0.7,0.3,0.4,0.4, 0.6),
                    nrow = 2, byrow = TRUE)



DatCases<- tibble(dgp_id= paste0("case", 1:4, ""), 
                   sid = I(list(c("A1,B1"), "B2", c("A2", "C1"), "C2")),
                   Ktrue = rep(2,4),  p = c(3,rep(5,3)), 
                   mu = I(list(mu_p3, mu_p5_all, mu_p5_three,mu_p5_three)), 
                   lambda1 = c(25,10, 25,10), 
                   lambda2 = c(35,5,35,5),
                   pi1 = rep(0.45, 4), 
                   pi2= rep(0.55, 4),
                   saliency = I(list(c(1,1,1),
                                     c(1,1,1,1,1),
                                     c(1,1,1,0,0),
                                     c(1,1,1,0,0))))

#Generate some plots
for (s in 1:4){
  K<- DatCases$Ktrue[s]
  p<-DatCases$p[s]
  mu_mat<-  DatCases$mu[[s]]
  lambda_vec <- c(DatCases$lambda1[s], DatCases$lambda2[s])
  p_vec<- c(DatCases$pi1[s], DatCases$pi2[s])
  dat<-generateData(n,p,K,mu_mat, lambda_vec, p_vec)
  x<-as.data.frame(dat[[1]])
  z<-dat[[2]]
  cluster<-factor(z)
  
  png(sprintf("simulations/figures/DGPs/caseGGpairs_%02d.png", s), width = 1200, height = 1000)
  p<-ggpairs(x, columns = 1:p,
             ggplot2::aes(color = cluster, alpha = 0.5),
             upper = list(continuous = wrap("cor", size = 3)),
             lower = list(continuous = wrap("points", alpha = 0.3))) +
    theme_bw(base_size = 10)
  
  print(p)
  dev.off()
}
#Generate 100 datasets for each scenario
R<-100 
n<- 2000
S<-4
datlistAll<-list()
alloclistAll<-list()

for (s in 1:S){
  datlist<-list()
  alloclist<-list()
  for (i in 1:R){
    ID<-DatCases$dgp_id[s]
    K<- DatCases$Ktrue[s]
    p<-DatCases$p[s]
    mu_mat<-  DatCases$mu[[s]]
    lambda_vec <- c(DatCases$lambda1[s], DatCases$lambda2[s])
    p_vec<- c(DatCases$pi1[s], DatCases$pi2[s])
    dat<-generateData(n,p,K,mu_mat, lambda_vec, p_vec)
    
    #Separate the data (x) and true allocation (z)
    datlist[[i]]<- dat[[1]]
    alloclist[[i]]<-dat[[2]]
  }
  datlistAll[[s]]<- datlist
  alloclistAll[[s]]<-alloclist
}
names(datlistAll)<- DatCases$dgp_id
names(alloclistAll)<- DatCases$dgp_id


# datlistALL:large list of 4 lists
#case1:List of 100 :generated data for K = 2, p=3 and n=2000
#case2:List of 100 :generated data for K = 2, p=5 and n=2000
#case3:List of 100 :generated data for K = 2, p=5 and n=2000
#case4:List of 100 :generated data for K = 2, p=5 and n=2000 (Feature Saliency)


#Save these into a .Rdata to be used in the simulations
saveRDS(datlistAll, file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_dpgs1_4.rds")
saveRDS(list(DatCases,alloclistAll), file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_dpgs1_4.rds")



####################################################################################
#Changes made on 4/3/26
####################################################################################
set.seed(2025)
#-----------DATA CASES-----------------
n<- 2000
#All variables informative
mu_p5_all<-matrix(c(0.2, 0.5, 0.8,0.4,0.5, 
                    0.7, 0.3, 0.4, 0.5, 0.9),
                  nrow = 2, byrow = TRUE)

mu_p10_all <- matrix(c(
  # Cluster 1 (10 variables)
  0.2, 0.5, 0.8, 0.3, 0.5, 0.3, 0.6, 0.7, 0.2, 0.9,
  # Cluster 2 (10 variables)
  0.7, 0.3, 0.4, 0.5, 0.9, 0.8, 0.2, 0.1, 0.6, 0.4
), nrow = 2, byrow = TRUE)

#5 variables non-informative
mu_p10_five<- matrix(c(
  0.2, 0.5, 0.8, 0.3, 0.5, 0.3, 0.6, 0.7, 0.2, 0.4,
  0.7, 0.3, 0.4, 0.5, 0.9, 0.3, 0.6, 0.7, 0.2, 0.4
), nrow = 2, byrow = TRUE)


DatCases<- tibble(dgp_id= paste0("case", 1:4, ""), 
                  sid = I(list(c("A1,B1"), "B2", c("A2", "C1"), "C2")),
                  Ktrue = rep(2,4),  p = c(5,rep(10,3)), 
                  mu = I(list(mu_p5_all,mu_p10_all,mu_p10_five,mu_p10_five)), 
                  lambda1 = c(25,10, 25,10), 
                  lambda2 = c(35,5,35,5),
                  pi1 = rep(0.45, 4), 
                  pi2= rep(0.55, 4),
                  saliency = I(list(c(1,1,1,1,1),
                                    rep(1,10),
                                    c(rep(1,5), rep(0,5)),
                                    c(rep(1,5), rep(0,5)))))

#Generate some plots
for (s in 1:4){
  K<- DatCases$Ktrue[s]
  p<-DatCases$p[s]
  mu_mat<-  DatCases$mu[[s]]
  lambda_vec <- c(DatCases$lambda1[s], DatCases$lambda2[s])
  p_vec<- c(DatCases$pi1[s], DatCases$pi2[s])
  dat<-generateData(n,p,K,mu_mat, lambda_vec, p_vec)
  x<-as.data.frame(dat[[1]])
  z<-dat[[2]]
  cluster<-factor(z)
  
  png(sprintf("simulations/figures/DGPs/caseGGpairs_higherDim_%02d.png", s), width = 1200, height = 1000)
  p<-ggpairs(x, columns = 1:p,
             ggplot2::aes(color = cluster, alpha = 0.5),
             upper = list(continuous = wrap("cor", size = 3)),
             lower = list(continuous = wrap("points", alpha = 0.3))) +
    theme_bw(base_size = 10)
  
  print(p)
  dev.off()
}
#Generate 100 datasets for each scenario
R<-100 
S<-4
datlistAll<-list()
alloclistAll<-list()

for (s in 1:S){
  datlist<-list()
  alloclist<-list()
  for (i in 1:R){
    ID<-DatCases$dgp_id[s]
    K<- DatCases$Ktrue[s]
    p<-DatCases$p[s]
    mu_mat<-  DatCases$mu[[s]]
    lambda_vec <- c(DatCases$lambda1[s], DatCases$lambda2[s])
    p_vec<- c(DatCases$pi1[s], DatCases$pi2[s])
    dat<-generateData(n,p,K,mu_mat, lambda_vec, p_vec)
    
    #Separate the data (x) and true allocation (z)
    datlist[[i]]<- dat[[1]]
    alloclist[[i]]<-dat[[2]]
  }
  datlistAll[[s]]<- datlist
  alloclistAll[[s]]<-alloclist
}
names(datlistAll)<- DatCases$dgp_id
names(alloclistAll)<- DatCases$dgp_id


# datlistALL:large list of 4 lists
#case1:List of 100 :generated data for K = 2, p=5 and n=2000
#case2:List of 100 :generated data for K = 2, p=10 and n=2000
#case3:List of 100 :generated data for K = 2, p=10 and n=2000
#case4:List of 100 :generated data for K = 2, p=10 and n=2000 (Feature Saliency)


#Save these into a .Rdata to be used in the simulations
saveRDS(datlistAll, file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_dpgs1_4_higherp.rds")
saveRDS(list(DatCases,alloclistAll), file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_dpgs1_4_higherp.rds")






















