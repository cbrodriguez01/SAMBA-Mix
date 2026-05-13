#===================================================================
# Simulate data cases
# Programmer: Carmen Rodriguez C.
# Date updated: 2025/12/15
#===================================================================
#setwd("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/")

library(dplyr)  # Data wrangling and manipulation
library(tidycensus)
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
mu_p5_all<-matrix(c(0.2, 0.5, 0.8,0.4,0.5, 
            0.7, 0.3, 0.4, 0.5, 0.9),
            nrow = 2, byrow = TRUE)

mu_p5_three<-matrix(c(0.2,0.3,0.8,0.4,  0.5, 
                      0.7,0.3,0.4,0.4, 0.9),
                  nrow = 2, byrow = TRUE)

DatCases<- tibble( cid= paste0("case", 1:5, ""),
  Ktrue = rep(2,5),  p = c(3,rep(5,4)), 
  overlap = c("none", "moderate", "moderate", "none", "moderate"), 
  mu = I(list(mu_p3, mu_p5_all, mu_p5_all,mu_p5_three,mu_p5_three)), 
  lambda1 = c(25,10,10, 25,10), 
  lambda2 = c(35,5,10,35,5),
  pi1 = rep(0.45, 5), 
  pi2= rep(0.55, 5),
  saliency = I(list(c(1,1,1),
                    c(1,1,1,1,1),
                    c(1,1,1,1,1),
                    c(1,0,1,0,1),
                    c(1,0,1,0,1))))

#Generate some plots
for (s in 1:5){
K<- DatCases$Ktrue[s]
p<-DatCases$p[s]
mu_mat<-  DatCases$mu[[s]]
lambda_vec <- c(DatCases$lambda1[s], DatCases$lambda2[s])
p_vec<- c(DatCases$pi1[s], DatCases$pi2[s])
dat<-generateData(n,p,K,mu_mat, lambda_vec, p_vec)
x<-as.data.frame(dat[[1]])
z<-dat[[2]]
cluster<-factor(z)

png(sprintf("simulations/figures/Cases1_5_Data/caseGGpairs_%02d.png", s), width = 1200, height = 1000)
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
S<-5
datlistAll<-list()
alloclistAll<-list()

for (s in 1:S){
  datlist<-list()
  alloclist<-list()
for (i in 1:R){
    ID<-DatCases$cid[s]
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
names(datlistAll)<- DatCases$cid
names(alloclistAll)<- DatCases$cid


# datlistALL:large list of 5 lists
    #case1:List of 100 :generated data for K = 2, p=3 and n=2000
    #case2:List of 100 :generated data for K = 2, p=5 and n=2000
    #case3:List of 100 :generated data for K = 2, p=5 and n=2000
    #case4:List of 100 :generated data for K = 2, p=5 and n=2000 (Feature Saliency)
    #case4:List of 100 :generated data for K = 2, p=5 and n=2000 (Feature Saliency)

#Save these into a .Rdata to be used in the simulations
saveRDS(datlistAll, file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_cases1_5.rds")
saveRDS(list(DatCases,alloclistAll), file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_cases1_5.rds")

#-----------MOCK ACS Data-----------------
#Project Brainstorming doc in Gdrive
#zero inflation and asymmetric (lower concentration parameters to exagerate skweness)
#Some variables at the boundary of (0,1)
#n was selected based on largest number of census tracts ( ~9000 for California), and census block group
#for some states if the desire is to  to explore how small we can get geographically
#We do the 100 samples of 2000, but do one large n exploration once we are settle on prior parameters
#simulate_ACS_data
#https://journal.r-project.org/articles/RJ-2015-019/
set.seed(1994)
n<-2000
Ktrue <-5
p<-20
pi_true1<-rep(0.2,5)
pi_true2<-c(0.02, rep(0.245,4)) #one very small component to test
lambda_vec<- runif(Ktrue, 3, 10)  #to have a range of skewed vs. concentrated
#Remember:  
#lower lambda: more skewed, diffuse
#larger lambda: tighter around mu

# All variables are important: Case 6
mu_mat_case6 <- matrix(NA, nrow = Ktrue, ncol = p) 

#Note these values were selected based on the distribution of ACS data we explored in a prev paper
for(k in 1:Ktrue){
  # right-skewed “ACS-like” vars -- lower "signal" like discussed with Briana on 12/15/25
  skew_vars <- c(1,2,4,5,6)
  mu_mat_case6[k, skew_vars] <- runif(length(skew_vars), 0.05, 0.15)
  
  
  mid_vars <- c(3,8,10,11)
  mu_mat_case6[k, mid_vars] <- runif(length(mid_vars),0.35,0.65)
  
  #high signals
  high_vars <- setdiff(1:p, c(skew_vars,mid_vars))
  mu_mat_case6[k, high_vars] <- runif(length(high_vars),0.70,0.995)
}

# 
# dat<-generateData(n,p,Ktrue,mu_mat_case6, lambda_vec,pi_true1)
# x<-as.data.frame(dat[[1]])
# z<-dat[[2]]
# cluster<-factor(z)
# 
# dat2<-generateData(n,p,Ktrue,mu_mat_case6, lambda_vec,pi_true2)
# x1<-as.data.frame(dat2[[1]])
# z1<-dat2[[2]]
# cluster1<-factor(z1)

# p1<-ggpairs(x, columns = c(1,2,4,5,6),
#         ggplot2::aes(color = cluster, alpha = 0.5),
#         upper = list(continuous = wrap("cor", size = 3)),
#         lower = list(continuous = wrap("points", alpha = 0.3))) +
#   theme_bw(base_size = 10)
# 
# 
# p2<-ggpairs(x, columns = c(3,8,10,11),
#         ggplot2::aes(color = cluster, alpha = 0.5),
#         upper = list(continuous = wrap("cor", size = 3)),
#         lower = list(continuous = wrap("points", alpha = 0.3))) +
#   theme_bw(base_size = 10)
# 
# 
# p3<-ggpairs(x, columns = c(7,9, 12, 13, 14, 15, 16, 17, 18, 19, 20),
#         ggplot2::aes(color = cluster, alpha = 0.5),
#         upper = list(continuous = wrap("cor", size = 3)),
#         lower = list(continuous = wrap("points", alpha = 0.3))) +
#   theme_bw(base_size = 10)
# 
# 
# p4<-ggpairs(x1, columns = c(1,2,4,5,6),
#         ggplot2::aes(color = cluster1, alpha = 0.5),
#         upper = list(continuous = wrap("cor", size = 3)),
#         lower = list(continuous = wrap("points", alpha = 0.3))) +
#   theme_bw(base_size = 10)
# 
# 
# p5<-ggpairs(x1, columns = c(3,8,10,11),
#         ggplot2::aes(color = cluster1, alpha = 0.5),
#         upper = list(continuous = wrap("cor", size = 3)),
#         lower = list(continuous = wrap("points", alpha = 0.3))) +
#   theme_bw(base_size = 10)
# 
# 
# p6<-ggpairs(x1, columns = c(7,9, 12, 13, 14, 15, 16, 17, 18, 19, 20),
#         ggplot2::aes(color = cluster1, alpha = 0.5),
#         upper = list(continuous = wrap("cor", size = 3)),
#         lower = list(continuous = wrap("points", alpha = 0.3))) +
#   theme_bw(base_size = 10)


##Export these for supplementary materials
# path<-"simulations/figures/CensusSim/"
# 
# pdf(file = paste0(path, "5clustACS_mock_same_pi.pdf"),   # The directory you want to save the file in
#     width = 10, # The width of the plot in inches
#     height = 20) # The height of the plot in inches
# p1
# p2
# p3
# dev.off()


# pdf(file = paste0(path, "5clustACS_mock_1small_pi.pdf"),   # The directory you want to save the file in
#     width = 10, # The width of the plot in inches
#     height = 20) # The height of the plot in inches
# p4
# p5
# p6
# dev.off()

# Only 12 variables are important: Case 7
#For the 8 unimportant variables we make the distribution identical across all components
#lambda_vec stays the same

#Randomly select which 8 variables are unimportant:
uindex<- sample(1:p,8)
mu_mat_case7 <- mu_mat_case6
global_mu_j<- 0.5
for (j in uindex){
  mu_mat_case7[, j] <- global_mu_j
}
  
# dat<-generateData(n,p,Ktrue,mu_mat_case7, lambda_vec,pi_true1)
# x<-as.data.frame(dat[[1]])
# z<-dat[[2]]
# cluster<-factor(z)
# 
# dat2<-generateData(n,p,Ktrue,mu_mat_case7, lambda_vec,pi_true2)
# x1<-as.data.frame(dat2[[1]])
# z1<-dat2[[2]]
# cluster1<-factor(z1)

# 
# p1<-ggpairs(x, columns = c(1,2,4,5,6),
#             ggplot2::aes(color = cluster, alpha = 0.5),
#             upper = list(continuous = wrap("cor", size = 3)),
#             lower = list(continuous = wrap("points", alpha = 0.3))) +
#   theme_bw(base_size = 10)
# 
# 
# p2<-ggpairs(x, columns = c(3,8,10,11),
#             ggplot2::aes(color = cluster, alpha = 0.5),
#             upper = list(continuous = wrap("cor", size = 3)),
#             lower = list(continuous = wrap("points", alpha = 0.3))) +
#   theme_bw(base_size = 10)
# 
# 
# p3<-ggpairs(x, columns = c(7,9, 12, 13, 14, 15, 16, 17, 18, 19, 20),
#             ggplot2::aes(color = cluster, alpha = 0.5),
#             upper = list(continuous = wrap("cor", size = 3)),
#             lower = list(continuous = wrap("points", alpha = 0.3))) +
#   theme_bw(base_size = 10)
# 
# 
# p4<-ggpairs(x1, columns = c(1,2,4,5,6),
#             ggplot2::aes(color = cluster1, alpha = 0.5),
#             upper = list(continuous = wrap("cor", size = 3)),
#             lower = list(continuous = wrap("points", alpha = 0.3))) +
#   theme_bw(base_size = 10)
# 
# 
# p5<-ggpairs(x1, columns = c(3,8,10,11),
#             ggplot2::aes(color = cluster1, alpha = 0.5),
#             upper = list(continuous = wrap("cor", size = 3)),
#             lower = list(continuous = wrap("points", alpha = 0.3))) +
#   theme_bw(base_size = 10)
# 
# 
# p6<-ggpairs(x1, columns = c(7,9, 12, 13, 14, 15, 16, 17, 18, 19, 20),
#             ggplot2::aes(color = cluster1, alpha = 0.5),
#             upper = list(continuous = wrap("cor", size = 3)),
#             lower = list(continuous = wrap("points", alpha = 0.3))) +
#   theme_bw(base_size = 10)
# 
# 
# ##Export these for supplementary materials
# path<-"simulations/figures/CensusSim/"
# pdf(file = paste0(path, "5clustACS_mock_same_pi_FS.pdf"),   # The directory you want to save the file in
#     width = 10, # The width of the plot in inches
#     height = 20) # The height of the plot in inches
# p1
# p2
# p3
# dev.off()
# pdf(file = paste0(path, "5clustACS_mock_1small_pi_FS.pdf"),   # The directory you want to save the file in
#     width = 10, # The width of the plot in inches
#     height = 20) # The height of the plot in inches
# p4
# p5
# p6
# dev.off()

#Generate 100 datasets for each scenario (case6 and case7)
#For the small pi, only do this for case6,no need to repeat for case7
datlistAll_6<-list()
alloclistAll_6<-list()
for (i in 1:R){
    dat<-generateData(n,p,Ktrue,mu_mat_case6, lambda_vec,pi_true1)
    #Separate the data (x) and true allocation (z)
    datlistAll_6[[i]]<- dat[[1]]
    alloclistAll_6[[i]]<-dat[[2]]}

datlistAll_7<-list()
alloclistAll_7<-list()
for (i in 1:R){
  dat<-generateData(n,p,Ktrue,mu_mat_case7, lambda_vec,pi_true1)
  #Separate the data (x) and true allocation (z)
  datlistAll_7[[i]]<- dat[[1]]
  alloclistAll_7[[i]]<-dat[[2]]}


dat_pitrue2<-generateData(n,p,Ktrue,mu_mat_case6, lambda_vec,pi_true2)
#Save true values
DatCases2<- tibble(cid= paste0("case", 6:7, ""),
                   Ktrue = rep(5,2),  p = c(20,20), 
                   overlap = c( "moderate",  "moderate"), 
                   mu = I(list(mu_mat_case6,mu_mat_case7)), 
                   lambda_vec = I(list(lambda_vec,lambda_vec)), 
                   pivec = I(list(pi_true1,pi_true1)), 
                   saliencyNum = c(20,12))


#Save these into a .Rdata to be used in the simulations
saveRDS(list(datlistAll_6,dat_pitrue2), file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_case6.rds")
saveRDS(list(DatCases2,alloclistAll_6), file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case6.rds")


saveRDS(datlistAll_7, file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_case7.rds")
saveRDS(alloclistAll_7, file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case7.rds")



