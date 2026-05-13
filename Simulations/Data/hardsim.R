library(tidyverse)
library(ggplot2)
library(GGally)
source("R/helpers.R")
source("R/bayesmbmm_OFMM_SF_v2.R")



#Stress test 
#K= 3
#p= 20
#n= 1500
# | Type             | Features | Description                       |
#   | ---------------- | -------- | --------------------------------- |
#   | Strong signal    | 1–6      | Clear separation                  |
#   | Weak signal      | 7–12     | Overlapping distributions         |
#   | Correlated noise | 13–16    | No cluster signal, but correlated |
#   | Pure noise       | 17–20    | Independent noise                 |

n <- 2000
K <- 3
p <- 20
pitrue<-c(0.4,0.35,0.25)
lambda_vec<-c(30,20,10) #make it difficult
mu <- matrix(NA, nrow = K, ncol = p)

#strong signal
mu[,1:6] <- rbind(
  rep(0.8, 6),
  rep(0.2, 6),
  rep(0.5, 6)
)

#weak signal
mu[,7:12] <- rbind(
  rep(0.55, 6),
  rep(0.45, 6),
  rep(0.50, 6)
)

#simulate “pseudo-correlation” via similar mu patterns:
base_patterns <- rbind(
  c(0.3, 0.5, 0.7),
  c(0.35, 0.55, 0.75),
  c(0.32, 0.52, 0.72),
  c(0.28, 0.48, 0.68)
)

for (j in 13:16) {
  mu[,j] <- base_patterns[j-12, ]
}

#pure noise
mu[,17:20] <- 0.5

dat<-generateData(n,p,K,mu, lambda_vec,pitrue)
x<-as.data.frame(dat[[1]])
z<-dat[[2]]
cluster<-factor(z)
ggpairs(x,ggplot2::aes(color = cluster, alpha = 0.5),
        upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 10)



# 
# fit <- Ombmm_FeatSaliency(
#     x = x,
#     K = 15,
#     niter = 10000,
#     nburnin = 5000,
#     nchain = 1,
#     init_list = NULL,
#     g = 5, h = 3,
#     a_mu0 = 1, b_mu0 = 1,
#     a_mu1 = 1, b_mu1 = 1,
#     a_rho = 5, b_rho = 1,
#     e0 = 0.1,
#     priorOnAlpha = "no",
#     plot_trace = TRUE
#   )

set.seed(1994)
R<-50
datlist<-list()
alloclist<-list()
for (i in 1:R){
  dat<-generateData(n,p,K,mu, lambda_vec,pitrue)
  #Separate the data (x) and true allocation (z)
  datlist[[i]]<- dat[[1]]
  alloclist[[i]]<-dat[[2]]}
out<-list(datlist, alloclist, K, mu, lambda_vec, pitrue)

saveRDS(out, file = "/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistStress.rds")


