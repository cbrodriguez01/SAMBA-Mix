# Simulate data mock acs data cases
# Programmer: Carmen Rodriguez C.
# Date updated: 2/16/2026
#===================================================================
#setwd("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/")

library(tidyverse)
library(ggplot2)
library(GGally)
source("R/helpers.R")

#-----------MOCK ACS Data-----------------
#Project Brainstorming doc in Gdrive
#zero inflation and asymmetric (lower concentration parameters to exagerate skweness)
#Some variables at the boundary of (0,1)
#for some states if the desire is to  to explore how small we can get geographically
#We do the 100 samples of 2000, but do one large n exploration once we are settle on prior parameters
#
#https://journal.r-project.org/articles/RJ-2015-019/
set.seed(2026)
n<-2000
Ktrue <-5
p<-12 #important vars
 #12 important + 6 noise vars
nk<-c(600, 500, 400, 300, 200)
pi_true1<- nk/n
pi_true2<-rep(0.2,5)
lambda_vec1 <- runif(5, 10,35) #[1] 19.87745 17.46102 10.38238 12.85730 17.44127
lambda_vec2<- rep(25,5)

# Staggered means (5 clusters x 12 variables)
# Dim 1-4: Income/Edu, Dim 5-8: Age/HH Size, Dim 9-12: Housing/Geography
mu_mat <- matrix(c(
  0.8, 0.7, 0.8, 0.9,  0.2, 0.1, 0.2, 0.1,  0.5, 0.5, 0.4, 0.3, # Cluster 1: High SES, Young
  0.2, 0.1, 0.2, 0.1,  0.8, 0.9, 0.7, 0.8,  0.1, 0.1, 0.2, 0.2, # Cluster 2: Low SES, Older
  0.5, 0.5, 0.6, 0.5,  0.5, 0.6, 0.5, 0.4,  0.8, 0.9, 0.7, 0.8, # Cluster 3: Middle, Suburban
  0.1, 0.2, 0.1, 0.2,  0.2, 0.1, 0.3, 0.2,  0.3, 0.2, 0.1, 0.2, # Cluster 4: Low SES, Young
  0.9, 0.8, 0.7, 0.8,  0.8, 0.7, 0.9, 0.8,  0.7, 0.6, 0.8, 0.9  # Cluster 5: High SES, High Age
), nrow = 5, byrow = TRUE)

#check distance to make sure there is not a lot of overlap
#recommended to add a small random perturbation to avoid deterministic separation
D <- as.matrix(dist(mu_mat))
round(D, 3)
colMeans(mu_mat)

##Generate data
# dat<-generateData(n,p,Ktrue,mu_mat, lambda_vec1,pi_true1)
# x<-as.data.frame(dat[[1]])
# z<-dat[[2]]
# cluster<-factor(z)
# 
# dat1<-generateData(n,p,Ktrue,mu_mat, lambda_vec2,pi_true2)
# x1<-as.data.frame(dat1[[1]])
# z1<-dat1[[2]]
# cluster1<-factor(z1)
# 
# dat2<-generateData(n,p,Ktrue,mu_mat, lambda_vec1,pi_true2)
# x2<-as.data.frame(dat2[[1]])
# z2<-dat2[[2]]
# cluster2<-factor(z2)

# lambda all the same, and pi varies
dat<-generateData(n,p,Ktrue,mu_mat, lambda_vec2,pi_true1)
x<-as.data.frame(dat[[1]])
z<-dat[[2]]
cluster<-factor(z)

dat1<-generateData(n,p,Ktrue,mu_mat, lambda_vec1,pi_true1)
x1<-as.data.frame(dat1[[1]])
z1<-dat1[[2]]
cluster1<-factor(z1)



# #Lambda varies across, pi varies across
ggpairs(x, columns = 1:12,
            ggplot2::aes(color = cluster, alpha = 0.5),
            upper = list(continuous = wrap("cor", size = 3)),
            lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 10)
# 
# # lambda same across, pi all the same--This one looks good
# ggpairs(x1, columns = 1:12,
#             ggplot2::aes(color = cluster1, alpha = 0.5),
#             upper = list(continuous = wrap("cor", size = 3)),
#             lower = list(continuous = wrap("points", alpha = 0.3))) +
#   theme_bw(base_size = 10)
# 
# 
# #Good but harder to identify: lambda varies, pi all the same
# ggpairs(x2, columns = 1:12,
#         ggplot2::aes(color = cluster2, alpha = 0.5),
#         upper = list(continuous = wrap("cor", size = 3)),
#         lower = list(continuous = wrap("points", alpha = 0.3))) +
#   theme_bw(base_size = 10)
# 
# 
# #This is a mix of x1 and x2; i think this could be easy: lambda all the same, and pi varies
# ggpairs(x3, columns = 1:12,
#         ggplot2::aes(color = cluster3, alpha = 0.5),
#         upper = list(continuous = wrap("cor", size = 3)),
#         lower = list(continuous = wrap("points", alpha = 0.3))) +
#   theme_bw(base_size = 10)
# 
# 
# dat_list_mkacs1<- list(x,x1,x2,x3)

###############################################################################
###############################################################################
##Now generate data for p=20
# Only 12 variables are important
#For the 8 unimportant variables we make the distribution identical across all components

p1<-20
mu_mat2 <- matrix(NA, Ktrue, p1)
mu_mat2[, 1:12] <- mu_mat

global_mu_j<- 0.5
for (j in 13:20){
  mu_mat2[, j] <- global_mu_j
}

##Generate data
dat<-generateData(n,p1,Ktrue,mu_mat2, lambda_vec1,pi_true1)
x<-as.data.frame(dat[[1]])
z<-dat[[2]]
cluster<-factor(z)

dat1<-generateData(n,p1,Ktrue,mu_mat2, lambda_vec2,pi_true2)
x1<-as.data.frame(dat1[[1]])
z1<-dat1[[2]]
cluster1<-factor(z1)

dat2<-generateData(n,p1,Ktrue,mu_mat2, lambda_vec1,pi_true2)
x2<-as.data.frame(dat2[[1]])
z2<-dat2[[2]]
cluster2<-factor(z2)

dat3<-generateData(n,p1,Ktrue,mu_mat2, lambda_vec2,pi_true1)
x3<-as.data.frame(dat3[[1]])
z3<-dat3[[2]]
cluster3<-factor(z3)


#Lambda varies across, pi varies across
ggpairs(x, columns = 1:20,
        ggplot2::aes(color = cluster, alpha = 0.5),
        upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 10)

# lambda same across, pi all the same--This one looks good
ggpairs(x1, columns = 1:20,
        ggplot2::aes(color = cluster1, alpha = 0.5),
        upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 10)


#Good but harder to identify: lambda varies, pi all the same
ggpairs(x2, columns = 1:20,
        ggplot2::aes(color = cluster2, alpha = 0.5),
        upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 10)


#This is a mix of x1 and x2; i think this could be easy: lambda all the same, and pi varies
ggpairs(x3, columns = 1:20,
        ggplot2::aes(color = cluster3, alpha = 0.5),
        upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 10)


dat_list_mkacs2<- list(x,x1,x2,x3)


saveRDS(dat_list_mkacs1, file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistmockACS1.rds")
saveRDS(dat_list_mkacs2, file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistmockACS2.rds")


#Quick test!
source("R/bayesmbmm_OFMM_SF_v2.R")

start_time <- Sys.time()
#changing prior for lambda because we have different values
fit1<-Ombmm_FeatSaliency(
  x = x,
  K = 20,
  niter = 10000,
  nburnin = 0,
  nchain = 1,
  init_list = NULL, 
  g =  5, h = 3, a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
  a_rho = 1, b_rho = 20, e0 = 0.01,
  priorOnAlpha = "no",
  p.true = pi_true1, plot_trace = TRUE)

end_time <- Sys.time()


fit2<-Ombmm_FeatSaliency(
  x = x1,
  K = 20,
  niter = 10000,
  nburnin = 3000,
  nchain = 1,
  init_list = NULL, 
  g =  5, h = 3, a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
  a_rho = 1, b_rho = 20, e0 = 0.01,
  priorOnAlpha = "no",
  p.true = pi_true1, plot_trace = TRUE)

end_time <- Sys.time()
runtime_sec <- as.numeric(difftime(end_time, start_time, units = "secs")) # ~43 minutes






#Generate 100 datasets for each scenario (case6 and case7)
#For the small pi, only do this for case6,no need to repeat for case7
R<-100
datlist1<-list()
alloclist1<-list()
for (i in 1:R){
  dat<-generateData(n,p,Ktrue,mu_mat, lambda_vec,pi_true1)
  #Separate the data (x) and true allocation (z)
  datlist1[[i]]<- dat[[1]]
  alloclist1[[i]]<-dat[[2]]}

datlist2<-list()
alloclist2<-list()
for (i in 1:R){
  dat<-generateData(n,p1,Ktrue,mu_mat2, lambda_vec,pi_true1)
  #Separate the data (x) and true allocation (z)
  datlist2[[i]]<- dat[[1]]
  alloclist2[[i]]<-dat[[2]]}



#Save true values
DatCases2<- tibble(cid= paste0("mockACS", 1:2, ""),
                   Ktrue = rep(5,2),  p = c(12,20), 
                   mu = I(list(mu_mat, mu_mat2)), 
                   lambda_vec = I(list(lambda_vec,lambda_vec)), 
                   pivec = I(list(pi_true1,pi_true1)), 
                   saliencyNum = c(20,12))


#Save these into a .Rdata to be used in the simulations
saveRDS(list(datlist1,datlist2), file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistmockACS.rds")
saveRDS(list(DatCases2,alloclist1, alloclist2), file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistmockACS.rds")






##PLOTS


p1<-ggpairs(x1, columns = 1:5,
            ggplot2::aes(color = cluster, alpha = 0.5),
            upper = list(continuous = wrap("cor", size = 3)),
            lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 10)

p2<-ggpairs(x1, columns = 6:12,
            ggplot2::aes(color = cluster, alpha = 0.5),
            upper = list(continuous = wrap("cor", size = 3)),
            lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 10)



pdf(file = "simulations/figures/MockACSSimData/mock_acs1_022326.pdf")  
p1
p2
dev.off()



p3<-ggpairs(x1, columns = 1:5,
            ggplot2::aes(color = cluster1, alpha = 0.5),
            upper = list(continuous = wrap("cor", size = 3)),
            lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 10)

p4<-ggpairs(x1, columns = 6:12,
            ggplot2::aes(color = cluster1, alpha = 0.5),
            upper = list(continuous = wrap("cor", size = 3)),
            lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 10)


p5<-ggpairs(x1, columns = c(2,5,7,15),
            ggplot2::aes(color = cluster1, alpha = 0.5),
            upper = list(continuous = wrap("cor", size = 3)),
            lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 10)


pdf(file = "simulations/figures/MockACSSimData/mock_acs2_022326.pdf")  
p3
p4
p5
dev.off()







########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################


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
lambda_vec<- runif(Ktrue, 8,35)  #to have a range of skewed vs. concentrated (3,10)
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
  mu_mat_case6[k, mid_vars] <- runif(length(mid_vars),0.30,0.65)
  
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



#Added 1/19/26: Changes to the mock ACS data: case 5: 12 variables all important (the 12 vars from case 7), case 7: 20 variables where 8 are noise
important_index <- setdiff(1:p, uindex)
length(important_index)

mu_mat_case5 <- mu_mat_case7[, important_index, drop = FALSE]



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
R<-100
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

datlistAll_5<-list()
alloclistAll_5<-list()
for (i in 1:R){
  dat<-generateData(n,p=12,Ktrue,mu_mat_case5, lambda_vec,pi_true1)
  #Separate the data (x) and true allocation (z)
  datlistAll_5[[i]]<- dat[[1]]
  alloclistAll_5[[i]]<-dat[[2]]}


dat_pitrue2<-generateData(n,p,Ktrue,mu_mat_case6, lambda_vec,pi_true2)
#Save true values
DatCases2<- tibble(cid= paste0("case", c(6,7,5), ""),
                   Ktrue = rep(5,3),  p = c(20,20, 12), 
                   overlap = c( "moderate",  "moderate", "moderate"), 
                   mu = I(list(mu_mat_case6,mu_mat_case7, mu_mat_case5)), 
                   lambda_vec = I(list(lambda_vec,lambda_vec,lambda_vec)), 
                   pivec = I(list(pi_true1,pi_true1, pi_true1)), 
                   saliencyNum = c(20,12, 12))


#Save these into a .Rdata to be used in the simulations
saveRDS(list(datlistAll_6,dat_pitrue2), file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_case6.rds")
saveRDS(list(DatCases2,alloclistAll_6), file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case6.rds")


saveRDS(datlistAll_7, file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_case7.rds")
saveRDS(alloclistAll_7, file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case7.rds")


saveRDS(datlistAll_5, file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_case5.rds")
saveRDS(alloclistAll_5, file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case5.rds")



