#Hyperparameter sensitivity exploration just for the first 50 datasets cases 3 and 4
#This document includes: lambda, rho, e0, Kmax and prioronalpha
library(purrr)
library(tidyverse)
library(latex2exp)
files_path<-"/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/"


########################################################################################
########################################################################################
#Is the model sensitive to the prior on rho?
#Case 3 includes p=5 and saliency indicator: (1,1,1,0,0)
files <- list.files(files_path, full.names = TRUE, pattern = "^sensitivityrho_hp_case3_")
results_list <- lapply(files, function(f) readRDS(f))

results_df <- map_dfr(results_list, function(x) {
  cbind(
    rid = x$rid,
    x$hyperparams,
    as.data.frame(x$results)
  )
})

head(results_df)

#Calculate the E(rho), and use that in my plots
results_df<- results_df %>% mutate(E_rho = (a/(a+b)))

#Summarize ARI 
results_df %>% group_by(E_rho) %>%
  summarise(mean_ARI = mean(ARI, na.rm = TRUE),sd_ARI   = sd(ARI,na.rm = TRUE), .groups = "drop")


#This looks good; 
results_df %>%
  group_by(E_rho) %>%
  summarise(
    across(starts_with("PIP"), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )


results_df %>%
  group_by(E_rho) %>%
  summarise(
    prob_correct_K = mean(K == 2, na.rm = TRUE),
    mean_K = mean(K, na.rm = TRUE),
    .groups = "drop"
  )

# 
# results_df %>%
#   filter(is.na(PIP1) | is.na(PIP2) | is.na(PIP3) |
#            is.na(PIP4) | is.na(PIP5)) %>%
#   select(rid, a, b, K, starts_with("PIP"))
# 
# 
# 

#Make plot showing ARI-- NOT VERY INFORMATIVE!
ggplot(results_df, aes(x = factor(E_rho), y = ARI, fill = factor(E_rho))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.1, size = 1) +
  labs(
    x = TeX(r"($E(rho$))"),
    y = "ARI",
    title = ""
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


results_df_long<-results_df %>% pivot_longer(cols = paste0("PIP",1:10), names_to = "Feature", values_to = "PIPest")

results_df_long %>% ggplot(aes(x= Feature, y= PIPest, fill = Feature)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(~ E_rho) +
  labs(
    x = "Feature",
    y = "Posterior Inclusion Probabilities (PIP)",
    title = ""
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none" )

results_df_long$Feature1 <- factor(
  results_df_long$Feature,
  levels = unique(results_df_long$Feature),
  labels = paste0("Var", seq_along(unique(results_df_long$Feature)))
)

#results_df_long %>% filter(Feature =="PIP2")


P1<-results_df_long %>%  ggplot(aes(x = factor(E_rho), y = PIPest)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.8) +
  facet_wrap(~ Feature1, nrow = 1) +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    title = "",
    x = expression(E(rho)),
    y = "PIP"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none"
  )

df_summary <- results_df_long %>%
  group_by(E_rho, Feature1) %>%
  summarise(
    mean_pip = mean(PIPest, na.rm = TRUE),
    lower = quantile(PIPest, 0.25,na.rm = TRUE),
    upper = quantile(PIPest, 0.75,na.rm = TRUE),
    .groups = "drop"
  )

P2<-ggplot(df_summary, aes(x = E_rho, y = mean_pip)) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~ Feature1, nrow = 2) +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    title = "",
    x = expression(E(rho)),
    y = "Mean PIP"
  ) +
  theme_bw(base_size = 14)


png(filename = "~/bayesmbmm/simulations/figures/PIP_boxplot_sensitivity_hp.png", width = 800, height = 200)
P1
dev.off()

png(filename = "~/bayesmbmm/simulations/figures/PIP_mean_sensitivity_hp.png", width = 800, height = 400)
P2
dev.off()




# #Quick checks for NAs-- 
# datlist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_dpgs1_4.rds")
# alloclist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_dpgs1_4.rds")
# scenariosdat<- alloclist[[1]]
# case3 = datlist$case3
# true_params3<- get_true_params(scenariosdat, "case3")
# #true_params4<- get_true_params(scenariosdat, "case4")
# true.pi<- true_params3$true_pi
# #true_mu<- true_params4$true_mu
# p<-true_params3$numfeat
# alloc_true<- alloclist[[2]][["case3"]]
# 
# fit13 <- Ombmm_FeatSaliency(
#   x = case3[[13]],
#   K = 20,
#   niter = 1500,
#   nburnin = 500,
#   nchain = 1,
#   init_list = NULL, 
#   g =  5, h = 3, a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
#   a_rho = 2.5, b_rho = 7.5, e0 = 0.1, priorOnAlpha = "no",
#   p.true = true.pi, plot_trace = TRUE)
# 
# Kmode<-fit$K
# alloc_est<-fit$modeloutput_perchain$clusterMembershipPerMethod[1,] 
# ARI<- mclust::adjustedRandIndex(alloc_est, alloc_true[[rid]])
# PIP = as.vector(fit$FeatureSal_out[[2]])
# 
# 
# fit14 <- Ombmm_FeatSaliency(
#   x = case3[[14]],
#   K = 20,
#   niter = 1500,
#   nburnin = 500,
#   nchain = 1,
#   init_list = NULL, 
#   g =  5, h = 3, a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
#   a_rho = 1, b_rho = 9, e0 = 0.1, priorOnAlpha = "no",
#   p.true = true.pi, plot_trace = TRUE)
# 
# fit17 <- Ombmm_FeatSaliency(
#   x = case3[[17]],
#   K = 20,
#   niter = 1500,
#   nburnin = 500,
#   nchain = 1,
#   init_list = NULL, 
#   g =  5, h = 3, a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
#   a_rho = 0.5, b_rho = 9.5, e0 = 0.1, priorOnAlpha = "no",
#   p.true = true.pi, plot_trace = TRUE)
# 
# 
# fit4 <- Ombmm_FeatSaliency(
#   x = case3[[4]],
#   K = 20,
#   niter = 1500,
#   nburnin = 500,
#   nchain = 1,
#   init_list = NULL, 
#   g =  5, h = 3, a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
#   a_rho = 0.5, b_rho = 9.5, e0 = 0.1, priorOnAlpha = "no",
#   p.true = true.pi, plot_trace = TRUE)
# 
# fit8 <- Ombmm_FeatSaliency(
#   x = case3[[8]],
#   K = 20,
#   niter = 1500,
#   nburnin = 500,
#   nchain = 1,
#   init_list = NULL, 
#   g =  5, h = 3, a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
#   a_rho = 1, b_rho = 9, e0 = 0.1, priorOnAlpha = "no",
#   p.true = true.pi, plot_trace = TRUE)
# 



########################################################################################
########################################################################################
########################################################################################
########################################################################################
#Sensitivity for the	Dirichlet prior hyperparameters and Kmax 
#e0_vec <- c(0.1, 0.01, 0.005)
# priorOnAlpha_vec<- c("gam_05_05", "gam_1_2")
# Kmax_vec<-c(10,30)
# > hp
# Kmax    e0 priorOnAlpha e0_str prior_str                                        pattern
#   1    10 0.100           no  1e-01        no     Sensitivitye0_case4_K10_e0_1e-01_no\\.rds$
#   2    30 0.100           no  1e-01        no     Sensitivitye0_case4_K30_e0_1e-01_no\\.rds$
#   3    10 0.010           no  1e-02        no     Sensitivitye0_case4_K10_e0_1e-02_no\\.rds$
#   4    30 0.010           no  1e-02        no     Sensitivitye0_case4_K30_e0_1e-02_no\\.rds$
#   5    10 0.005           no  5e-03        no     Sensitivitye0_case4_K10_e0_5e-03_no\\.rds$
#   6    30 0.005           no  5e-03        no     Sensitivitye0_case4_K30_e0_5e-03_no\\.rds$
#   7    10    NA    gam_05_05     NA gam_05_05 Sensitivitye0_case4_K10_e0_NA_gam_05_05\\.rds$
#   8    30    NA    gam_05_05     NA gam_05_05 Sensitivitye0_case4_K30_e0_NA_gam_05_05\\.rds$
#   9    10    NA      gam_1_2     NA   gam_1_2   Sensitivitye0_case4_K10_e0_NA_gam_1_2\\.rds$
#   10   30    NA      gam_1_2     NA   gam_1_2   Sensitivitye0_case4_K30_e0_NA_gam_1_2\\.rds$
# 

files1 <- list.files(
  files_path,
  full.names = TRUE,
  pattern = "^Sensitivitye0_hp_case4_K[0-40]+_e0_.*_rid_[0-9]+\\.rds$"
)

results_list <- lapply(files1, function(f) readRDS(f))

results_to_df <- function(results_list) {
  bind_rows(lapply(results_list, function(x) {
    n <- length(x$Kmode)
    
    tibble(
      rep = seq_len(n),
      K_hat = x$Kmode,
      MSE_mu = x$MSE_mu,
      ARI = x$ARIest,
      Kmax = x$hyperparams$Kmax,
      e0 = x$hyperparams$e0,
      priorOnAlpha = x$hyperparams$priorOnAlpha
    )
  }))
}

df_reps <- results_to_df(results_list)

#Summarize
df_reps %>%
  group_by(Kmax, e0, priorOnAlpha) %>%
  summarise(
    K_mean = mean(K_hat, na.rm = TRUE),
    K_sd = sd(K_hat, na.rm = TRUE),
    ARI_mean = mean(ARI, na.rm = TRUE)
  )


df_reps1 <- df_reps  %>%
  mutate(
    e0_label = ifelse(is.na(e0), priorOnAlpha, paste0("e0=", e0)))
   # combo = paste0("Kmax=", Kmax, ", ", e0_label))


P1<-ggplot(df_reps1, aes(x = e0_label, y = K_hat)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
  facet_wrap(~ Kmax, scales = "free_x", labeller= labeller(Kmax = c( "10" = "Kmax = 10", "30" = "Kmax = 30"))) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  labs(
    x = expression("Hyperparameter setting " ~ e[0]),
    y = expression(hat(K)[0])) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




P2<-ggplot(df_reps1, aes(x = e0_label, y = MSE_mu)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
  facet_wrap(~ Kmax, scales = "free_x", labeller= labeller(Kmax = c( "10" = "Kmax = 10", "30" = "Kmax = 30"))) +
  labs(
    x = "Hyperparameter setting",
    y = "MSE for cluster-means") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


png(filename = "~/bayesmbmm/simulations/figures/Khat_Ke0alpha_sensitivity_hp.png", width = 800, height = 300)
P1
dev.off()

png(filename = "~/bayesmbmm/simulations/figures/MSEmu_Ke0alpha_sensitivity_hp.png", width = 800, height = 300)
P2
dev.off()



saveRDS(list(df_reps1, results_df), file = "/n/home03/crodriguezcabrera/Mbeta_Project/Sim_results/hyperparam_sens_e0_rho.rds")


########################################################################################
########################################################################################
#Sensitivity for the	concentration parameter prior--> ~ 2 hours

files2 <- list.files(
  files_path,
  full.names = TRUE,
  pattern = "^Sensitivitylambda_hp_.*\\.rds$"
)



test<-readRDS(files2[1])

results_list <- lapply(files2, function(f) {
  out <- readRDS(f)
  
  tibble(
    scen = out$hyperparams$scen,
    rep  = out$hyperparams$rep,
    g    = out$hyperparams$g,
    h    = out$hyperparams$h,
    Kmode = out$Kmode,
    lambda_est = I(list(out$lambda_est)),
    mse_lambda = out$mse_lambda,
    mab_lambda = out$mab_lambda,
    mse_mu = out$MSE_mu,
    ARI = out$ARI
  )
})

results_df <- bind_rows(results_list)

#str(results_df)
#, na.rm = TRUE
summary_df <- results_df %>%
  group_by(scen, g, h) %>%
  summarise(
    mean_mse_lambda = mean(mse_lambda, na.rm = TRUE),
    sd_mse_lambda   = sd(mse_lambda, na.rm = TRUE),
    mean_mab_lambda = mean(mab_lambda, na.rm = TRUE),
    sd_mab_lambda   = sd(mab_lambda, na.rm = TRUE),
    mean_K   = mean(Kmode, na.rm = TRUE),
    sd_K     = sd(Kmode, na.rm = TRUE),
    mean_ARI = mean(ARI, na.rm = TRUE),
    mean_mse_mu = mean(mse_mu, na.rm = TRUE),
    sd_mean_mu = sd(mse_mu, na.rm = TRUE),
    .groups = "drop"
  )
summary_df_other<- results_df %>% filter(is.na(mse_lambda))
#15/700; 1 with kmode = 3 and 14 with kmode = 1


#Checking the error files:Warning message:
#In Ombmm_FeatSaliency(x = casedat[[rid]], K = 10, niter = 3000,  :
   #                     Degenerate fit: Kmode < 2. Returning NULL.
                      






