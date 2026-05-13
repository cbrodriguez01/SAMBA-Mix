library(BayesBinMix)
library(coda)


#alloclist5<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case5.rds")
#alloclist7<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case7.rds")
#alloclist6<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case6.rds")

#scenariosdat<- alloclist6[[1]] #tibble of true parameters (mu, lambda, pi)
alloclist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/MockACSUpdated/alloclistAll_mockacsdpgs1_5.rds")
true_allocation<-alloclist[[2]]




files_path<-"/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/model_comp/"
files1 <- list.files(files_path, full.names = TRUE, pattern = "^BayesBinMix_ACSMOCK_MockACS3.*\\.rds$")
files2 <- list.files(files_path, full.names = TRUE, pattern = "^BayesBinMix_ACSMOCK_MockACS3.*\\.rds$")

results_list_case5 <- lapply(files1, function(f) readRDS(f)) #mm this is 82, need to check which reps failed and re-submit them
results_list_case7 <- lapply(files2, function(f) readRDS(f)) #mm this is 82, need to check which reps failed and re-submit them

#3/21/26: 77


get_mode <- function(x) {
  uniq_x <- unique(x)
  uniq_x[which.max(tabulate(match(x, uniq_x)))]
}

summarize_bayesbinmix_output<-function(fits, true_alloc){
  map_dfr(seq_along(fits), function(r) {
    fit <- fits[[r]]
    z_true <- true_alloc[[r]]
    z_hat  <- fit$clusterMembershipPerMethod[,2]
    K_hat<-get_mode(fit$K.mcmc)
    
    tibble(rep = r,
           ARI = adjustedRandIndex(z_true, z_hat),
           K_hat = K_hat,
           AcceptanceR= fit$chainInfo[2])
  })
}
res_case5 <- summarize_bayesbinmix_output(
  fits = results_list_case5,
  true_alloc = true_allocation$MockACS3)

res_case7 <- summarize_bayesbinmix_output(
  fits = results_list_case7,
  true_alloc = true_allocation$MockACS5)

res_case5$scen<- "MockACS3"
res_case7$scen<-"MockACS5"

res_all <- bind_rows(res_case5, res_case7)
res_all$K_true<-5
res_all$method<- "BayesBinMix"

#Average
res_all %>%
  group_by(as.factor(scen)) %>%
  summarise(
    mean_ARI = mean(ARI),
    sd_ARI   = sd(ARI),
    mean_K   = mean(K_hat),
    P_correct_K = mean(K_hat == K_true),
    P_over  = mean(K_hat > K_true),
    P_under = mean(K_hat < K_true))


ggplot(res_all, aes(x = ARI, fill = as.factor(scen))) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  facet_wrap(~ as.factor(scen)) +
  theme_minimal()

ggplot(res_all, aes(x = K_hat)) +
  geom_bar() +
  facet_wrap(~ as.factor(scen)) +
  geom_vline(
    aes(xintercept = K_true),
    linetype = "dashed",
    color = "red"
  ) +
  theme_minimal()



ggplot(res_all, aes(x = K_hat, y = ARI)) +
  geom_jitter(
    width = 0.25,
    height = 0,
    alpha = 0.7
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ as.factor(scen)) +
  theme_minimal(base_size = 13) +
  labs(
    x = expression(hat(K)),
    y = "Adjusted Rand Index (ARI)",
    title = "BayesBinMix: ARI vs Estimated Number of Clusters",
    subtitle = "median-binarized data"
  )


#Tries two delta temps to assure good mixing by looking at the acceptance rates or swap moves
# Mean ARI is about <0.001,mean K_hat is about 11-12, 0% correct recovery of K and K_hat is overestimated
#identical behavior for both temps--which yielded good acceptance rates and thus we can attest good mixing
#   -The poor performance is due to model misspecification induced by binarization, not due to poor mixing or inadequate tempering.
# The median threshold erases: skewness, tail behavior
# Explains heterogeneity via mode clusters-- many small clusters
#BayesBinMix treats variables as conditionally independent Bernoullis--The model cannot recover structure that no longer exists in the data.



#Figure out which reps failed and run again


present1 <- as.integer(
  sub(".*rep([0-9]+)\\.rds", "\\1", files1)
)

present2 <- as.integer(
  sub(".*rep([0-9]+)\\.rds", "\\1", files2)
)

missing1 <- setdiff(1:100, present1)
missing2 <- setdiff(1:100, present2)
missing1
missing2
#okay they're the same


#Check the data
#sim_data_binary<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/sim_data_binary_case5_case7.rds")

# Example: check rep 24
#dat <- sim_data_binary$case5[[24]]
#apply(dat, 2, function(x) length(unique(x)))



#Save results to home folder
saveRDS(res_all,file = "/n/home03/crodriguezcabrera/Mbeta_Project/Sim_results/bayesbinmix_sim_new.rds")





  