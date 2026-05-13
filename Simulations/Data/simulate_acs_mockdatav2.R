# Simulate data mock acs data cases again porque la otra version me tiene jartaaa
# Programmer: Carmen Rodriguez C.
# Date updated: 3/10/26
#===================================================================
library(tidyverse)
library(ggplot2)
library(GGally)
source("R/helpers.R")
#source("R/bayesmbmm_OFMM_SF_v2.R")

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
p<-12
pi_true<-rep(1/Ktrue,Ktrue)
#Remember:  
#lower lambda: more skewed, diffuse (more variability/spread)
#larger lambda: tighter around mu

# All variables are important
mu_mat <- matrix(runif(60, 0.05, 0.45), Ktrue, p)
#combination of high, medium and low signal
mu_mat[1, 1:3]   <- 0.85
mu_mat[2, 4:6]   <- 0.85
mu_mat[3, 7:9]   <- 0.85
mu_mat[4, 10:12] <- 0.85
mu_mat[5, c(1,5,9)] <- 0.85


#Dataset 1: fix lambda at 60
lambda_vec<- rep(60,Ktrue)

dat<-generateData(n,p,Ktrue,mu_mat, lambda_vec,pi_true)
x<-as.data.frame(dat[[1]])
z<-dat[[2]]
cluster<-factor(z)

#columns = c(1,2,3)
p1<-ggpairs(x,ggplot2::aes(color = cluster, alpha = 0.5),
            upper = list(continuous = wrap("cor", size = 3)),
            lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 10)

#Dataset 2: fix lambda at 35
lambda_vec1<- rep(35,Ktrue)
dat1<-generateData(n,p,Ktrue,mu_mat, lambda_vec1,pi_true)
x1<-as.data.frame(dat1[[1]])
z1<-dat1[[2]]
cluster1<-factor(z1)



p2<-ggpairs(x1, 
        ggplot2::aes(color = cluster1, alpha = 0.5),
        upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 10)



#Dataset 3: lambda varies
lambda_vec2<-c(10,18,25,35,60) 
dat2<-generateData(n,p,Ktrue,mu_mat, lambda_vec2,pi_true)
x2<-as.data.frame(dat2[[1]])
z2<-dat2[[2]]
cluster2<-factor(z2)

p3<-ggpairs(x2,
        ggplot2::aes(color = cluster2, alpha = 0.5),
        upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 10)




##Export these for supplementary materials
path<-"simulations/figures/MockACSSimData/"

pdf(file = paste0(path, "5clustACS_mock_same_pi_3.10.26.pdf"),   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 20) # The height of the plot in inches
p1
p2
p3
dev.off()


# Only 12 variables are important: Case 7
#For the 8 unimportant variables we make the distribution identical across all components
#lambda_vec stays the same
# noise part (8 variables)
mu_noise <- matrix(0.5, Ktrue, 8)
# combine
mu_mat1 <- cbind(mu_mat, mu_noise)

dat3<-generateData(n,p=20,Ktrue,mu_mat1, lambda_vec,pi_true)
x3<-as.data.frame(dat3[[1]])
z3<-dat3[[2]]
cluster3<-factor(z3)

dat4<-generateData(n,p=20,Ktrue,mu_mat1, lambda_vec2,pi_true)
x4<-as.data.frame(dat4[[1]])
z4<-dat2[[2]]
cluster4<-factor(z4)


p4<-ggpairs(x3, 
            ggplot2::aes(color = cluster3, alpha = 0.5),
            upper = list(continuous = wrap("cor", size = 3)),
            lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 10)


p5<-ggpairs(x4, 
            ggplot2::aes(color = cluster4, alpha = 0.5),
            upper = list(continuous = wrap("cor", size = 3)),
            lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 10)


pdf(file = paste0(path, "5clustACS_mock_pi_FS_3.10.26.pdf"),   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 20) # The height of the plot in inches
p4
p5

dev.off()


#Save true values
DatCases<- tibble(dgp_id= paste0("MockACS", 1:5, ""),
                   Ktrue = rep(5,5),  p = c(rep(12,3), 20,20), 
                   overlap = c( "low",  "moderate", "high", "low", "high"), 
                   mu = I(list(mu_mat,mu_mat, mu_mat, mu_mat1, mu_mat1)), 
                   lambda_vec = I(list(lambda_vec,lambda_vec1,lambda_vec2,lambda_vec, lambda_vec2)), 
                   pivec = I(list(pi_true,pi_true, pi_true,pi_true,pi_true)))


#Generate 100 datasets for each scenario (case6 and case7)
#For the small pi, only do this for case6,no need to repeat for case7
R<-100
S<-5
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
    lambda_vec <- DatCases$lambda_vec[[s]]
    p_vec<- DatCases$pivec[[s]]
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







# 
# #Large lambda = 60
# datlistAll_macs1<-list()
# alloclistAll_macs1<-list()
# for (i in 1:R){
#   #large lambda
#   dat<-generateData(n,p,Ktrue,mu_mat, lambda_vec,pi_true)
#   #Separate the data (x) and true allocation (z)
#   datlistAll_macs1[[i]]<- dat[[1]]
#   alloclistAll_macs1[[i]]<-dat[[2]]}
# 
# #Moderate lambda = 35
# datlistAll_macs2<-list()
# alloclistAll_macs2<-list()
# for (i in 1:R){
#   dat<-generateData(n,p,Ktrue,mu_mat, lambda_vec1,pi_true)
#   #Separate the data (x) and true allocation (z)
#   datlistAll_macs2[[i]]<- dat[[1]]
#   alloclistAll_macs2[[i]]<-dat[[2]]}
# 
# #Varying lambda = c(10,18,25,35,60)
# datlistAll_macs3<-list()
# alloclistAll_macs3<-list()
# for (i in 1:R){
#   dat<-generateData(n,p,Ktrue,mu_mat, lambda_vec2,pi_true)
#   #Separate the data (x) and true allocation (z)
#   datlistAll_macs3[[i]]<- dat[[1]]
#   alloclistAll_macs3[[i]]<-dat[[2]]}
# 
# ##---With noise vars---
# #Large lambda = 60
# datlistAll_macs4<-list()
# alloclistAll_macs4<-list()
# for (i in 1:R){
#   #large lambda
#   dat<-generateData(n,p=20,Ktrue,mu_mat1, lambda_vec,pi_true)
#   #Separate the data (x) and true allocation (z)
#   datlistAll_macs4[[i]]<- dat[[1]]
#   alloclistAll_macs4[[i]]<-dat[[2]]}
# 
# #Varying lambda = c(10,18,25,35,60)
# datlistAll_macs5<-list()
# alloclistAll_macs5<-list()
# for (i in 1:R){
#   dat<-generateData(n,p=20,Ktrue,mu_mat1, lambda_vec2,pi_true)
#   #Separate the data (x) and true allocation (z)
#   datlistAll_macs5[[i]]<- dat[[1]]
#   alloclistAll_macs5[[i]]<-dat[[2]]}



#Save these into a .Rdata to be used in the simulations
saveRDS(datlistAll, file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/MockACSUpdated/datlistAll_mockacsdpgs1_5.rds")
saveRDS(list(DatCases,alloclistAll), file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/MockACSUpdated/alloclistAll_mockacsdpgs1_5.rds")





