#Revised 1/22/26
#mclust output explore
library(mclust)
library(dplyr)
library(purrr)
library(ggplot2)
#output_mclust<-readRDS("/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/mclust_fits_cases5and7.rds")
output_mclust1<-readRDS("/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/mclust_fits_mockacs_new.rds")

datlist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/MockACSUpdated/datlistAll_mockacsdpgs1_5.rds") 
alloclist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/MockACSUpdated/alloclistAll_mockacsdpgs1_5.rds")

# alloclist5<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case5.rds")
# alloclist7<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case7.rds")
# alloclist6<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case6.rds")

acsmock3<-datlist$MockACS3
acsmock5<-datlist$MockACS5



scenariosdat<- alloclist[[1]] #tibble of true parameters (mu, lambda, pi)
true_allocation<-alloclist[[2]]

#case7_alloc<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case7.rds")

mclust_fits_mock3<-output_mclust1$mockacs3
mclust_fits_mock5<-output_mclust1$mockacs5


summarize_mclust_output<-function(fits, true_alloc,sid){
  map_dfr(seq_along(fits), function(r){
    fit <- fits[[r]]
    z_true <- true_alloc[[r]]
    z_hat  <- fit$classification
    
    tibble(case = sid,
            rep = r,
            ARI = adjustedRandIndex(z_true, z_hat),
            K_hat = fit$G_hat,
            model= fit$modelName)
  })
}

  
res_case3 <- summarize_mclust_output(
  fits = mclust_fits_mock3,
  true_alloc = true_allocation$MockACS3,
   sid = "MockACS3"
)

res_case5 <- summarize_mclust_output(
  fits = mclust_fits_mock5,
  true_alloc = true_allocation$MockACS5,
  sid = "MockACS5"
)
  
res_all <- bind_rows(res_case3, res_case5)
res_all$K_true<-5
res_all$method<- "mclust"

res_all %>%
  group_by(case) %>%
  summarise(
    mean_ARI = mean(ARI),
    sd_ARI   = sd(ARI),
    mean_K   = mean(K_hat),
    P_correct_K = mean(K_hat == K_true),
    P_over  = mean(K_hat > K_true),
    P_under = mean(K_hat < K_true))


ggplot(res_all, aes(x = ARI, fill = case)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  facet_wrap(~ case) +
  theme_minimal()



ggplot(res_all, aes(x = K_hat)) +
  geom_bar() +
  facet_wrap(~ case) +
  geom_vline(
    aes(xintercept = K_true),
    linetype = "dashed",
    color = "red"
  ) +
  theme_minimal()

#Summary
#1.When all variables carry signal, the logit-Gaussian approximation is tttt“good enough” to recover clusters.
# this is an optimistic setting -- clusters are well separated and so the logit transformaiton does not severely
#distort geometry.-- that is there is no irrelevant noise to distort covariance estimation
#

#2. K_hat is biased upward (mean about 7)
# ARI has large variability
# Creates more clusters to explain the variance (Spurious covariance structures --> extra components)
#variable importance is not orthogonal after logit... that is variables that are uninformative on the beta 
#scale are not uninformative on the logit scale.


#Save results to home folder
saveRDS(res_all,file = "/n/home03/crodriguezcabrera/Mbeta_Project/Sim_results/mclust_sim_new.rds")


