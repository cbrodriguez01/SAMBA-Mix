#=================================================================
# Used simulated data from Cases 1-4 (see excel file with values)
#For each case,  we generated 100 datasets and calculated MSE, MAB and ARI for all model parameters
#In this script, we combine all outputs after cluster run and summarize results
# Programmer: Carmen Rodriguez C.
# Date updated: 2/4/2026
#===========================================================
library(ggplot2)
library(ggridges)
library(tidyverse)
library(purrr)
library(patchwork)
library(hrbrthemes)

#------First we look at model parameters for all
#Okay there were 13 cases where Kmode = 1, and hence the model stopped because K > 1;
#I have to check if this happened for models where induced cluster overlaps

files_path<-"/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs"
#files <- list.files(files_path, full.names = TRUE, pattern = "^Ombmmresult_case[1-4].*\\.rds$")
#Quick checks
#results_list_Kt[192] #kmode is 3, replicate 91, case2
#results_list_Kt[179] #kmode 3, replicate 8, case2
#results_list_Kt <- results_list_Kt[-c(179,192)]
#new models with higher dimensions--- added 4/3/26
files <- list.files(files_path, full.names = TRUE, pattern = "^Ombmmresult_hp_case[1-4].*\\.rds$")


results_list_Kt <- lapply(files, function(f) readRDS(f)[[1]])#[[1]]Only keep the data.frames

#[[1]]Only keep the data.frames: these would be cases where

#lapply(results_list_Kt, function(df) sapply(df, class))
lapply(results_list_Kt, function(df) names(df)[sapply(df, is.character)])

results_list_Kt <- results_list_Kt[sapply(results_list_Kt, is.data.frame)]
results_all <- bind_rows(results_list_Kt)
str(results_all)
table(results_all$scenario) 
# case1 case2 case3 case4 
# 98    86    96   100 

summary_df <- results_all %>%
  group_by(as.factor(scenario)) %>%
  summarise(across(c(mse_pi, mse_mu, mse_lambda, mab_pi, mab_mu, mab_lambda, ARIval),
                   list(mean = mean, sd = sd), .names = "{.col}_{.fn}"))

#Save results to csv file
write_csv(summary_df, file = "~/Mbeta_Project/Sim_results/ombmm_dgp_1_4_hp.csv")




# Scenarios 1-4
results_all<- results_all %>% mutate(scenario_num = case_when(
  scenario == "case1" ~ 1,
  scenario == "case2" ~ 2,
  scenario == "case3" ~ 3,
  scenario == "case4" ~ 4
))

#Plots
mse_long <- pivot_longer(results_all, cols = starts_with("mse_"),
                         names_to = "parameter", values_to = "MSE")

mab_long <- pivot_longer(results_all, cols = starts_with("mab_"),
                         names_to = "parameter", values_to = "MAB")

param_labels <- c(
  mse_pi     = "pi",
  mse_lambda = "lambda",
  mse_mu     = "mu"
)
p1<-ggplot(mse_long, aes(x = factor(scenario_num), y = MSE, fill = factor(scenario_num))) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_wrap( ~ parameter,
    scales = "free_y",
    labeller = labeller(parameter = param_labels, .default = label_parsed)) +
  labs(title = "",
    x = "Scenario",
    y = "MSE"
  ) +
  theme_bw() +
  theme(strip.text = element_text(size = 16),legend.position = "none")


# p2<-ggplot(mse_long, aes(x = scenario, y = MSE, fill = scenario)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(width = 0.15, alpha = 0.25, size = 1, color = "black") +
#   scale_y_log10() +   # makes π panel readable
#   facet_wrap(
#     ~ parameter,
#     scales = "free_y",
#     labeller = labeller(parameter = param_labels, .default = label_parsed)
#   ) +
#   labs(
#     title = "MSE by data scenario for each parameter",
#     x = "Scenario",
#     y = "MSE (log scale)"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     strip.text = element_text(size = 14)
#   )

path<-"~/bayesmbmm/simulations/figures/"
jpeg(file = paste0(path, "OmbmmMSE_scen1_4_hp.jpeg"),  
  width = 800, 
  height = 400) 
p1
dev.off()

#Mean absolute bias
param_labels1 <- c(
  mab_pi     = "pi",
  mab_lambda = "lambda",
  mab_mu     = "mu"
)


p3<-ggplot(mab_long, aes(x = scenario, y = MAB, fill = scenario)) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_wrap( ~ parameter,
              scales = "free_y",
              labeller = labeller(parameter = param_labels1, .default = label_parsed)) +
  labs(title = "",
       x = "Scenario",
       y = "MAB"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 14), legend.position = "none"
  )



path<-"~/bayesmbmm/simulations/figures/"
jpeg(file = paste0(path, "OmbmmMAB_cases1_4_hp.jpeg"),  
     width = 800, 
     height = 400) 
p3
dev.off()


#Combine plots
mab_val<- mab_long %>% select(scenario, replicate,parameter, MAB) %>% rename(value = MAB)
mab_val$metric<- "MAB"
mab_val<- mab_val %>% mutate(params = sub("^(mab)_", "", parameter))

mse_val<- mse_long %>% select(scenario, replicate,parameter, MSE) %>% rename(value = MSE)
mse_val$metric<- "MSE"
mse_val<- mse_val %>% mutate(params = sub("^(mse)_", "", parameter))

combined<- rbind(mab_val, mse_val)


path<-"~/bayesmbmm/simulations/figures/"
jpeg(file = paste0(path, "OmbmmMABMSE_cases1_4_hp.jpeg"),  
     width = 800, 
     height = 400) 

ggplot(combined, aes(x = scenario, y = value, fill = scenario)) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_grid(metric ~ params,
              scales = "free_y",
              labeller = label_parsed) +
  scale_y_log10() +
  labs(title = "",
       x = "Scenario",
       y = "Value (log-scale)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 14), legend.position = "none"
  )

dev.off()




########################################################################################
########################################################################################
#---True and mean of the posterior distribution for each parameter
alloclist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_dpgs1_4_higherp.rds")
true_params<-alloclist[[1]]
str(true_params)
#estimated params
est_info <- results_all %>% select(scenario, replicate)
#nrow(est_info) == length(results_list_est1)

est_df <- lapply(files, function(f) {
  x <- readRDS(f)
  
  # skip bad runs
  if (is.null(x$est_params)) return(NULL)
  
  tibble(
    scenario = x[[1]]$scenario,
    replicate = x[[1]]$replicate,
    est = list(x$est_params)
  )
}) %>% bind_rows()


est_df_matched <- est_info %>%
  left_join(est_df, by = c("scenario", "replicate"))



#Get a tibble for each parameter
extract_estimates <- function(est_list, scenario, replicate) {
  pi_df <- tibble(
    scenario = scenario,
    replicate = replicate,
    parameter = paste0("pi", 1:2),
    estimate = est_list[[1]])
  
  lambda_df <- tibble(
    scenario = scenario,
    replicate = replicate,
    parameter = paste0("lambda", 1:2),
    estimate = est_list[[2]])
  
  # mu (K × p matrix)
  mu_mat <- est_list[[3]]
  K <- nrow(mu_mat)
  p <- ncol(mu_mat)
  
  mu_df <- expand.grid(
    component = 1:K,
    dim = 1:p) %>% mutate(parameter = paste0("mu", component, "_", dim),
      estimate = as.vector(mu_mat)
    ) %>% mutate(scenario = scenario,replicate = replicate) %>%
    select(scenario, replicate, parameter, estimate)
  
  bind_rows(pi_df, lambda_df, mu_df)
}


est_long <- purrr::pmap_dfr(
  est_df_matched,
  function(scenario, replicate, est) {
    if (is.null(est)) return(NULL)
    
    extract_estimates(
      est_list = est,
      scenario = scenario,
      replicate = replicate
    )
  }
)

#to match param names from truth and simulations
true_simple <- true_params %>%
  select(dgp_id, pi1, pi2, lambda1, lambda2) %>%
  pivot_longer(
    cols = pi1:lambda2,
    names_to = "parameter",
    values_to = "true_value"
  ) %>%
  rename(scenario = dgp_id)

#quick check
#as.vector(true_params$mu[[1]])

true_mu <- true_params %>%
  select(dgp_id, mu) %>%
  mutate(mu = map(mu, ~ {
    m <- .x
    expand.grid(component = 1:nrow(m), dim = 1:ncol(m)) %>%
      mutate(parameter = paste0("mu", component, "_", dim),
             true_value = as.vector(m))})) %>%
  unnest(mu) %>%
  rename(scenario = dgp_id) %>% select(scenario, parameter, true_value)

true_long <- bind_rows(true_simple, true_mu)
plot_df <- est_long %>%
  left_join(true_long, by = c("scenario", "parameter"))

str(plot_df)

# Scenarios 1-4
plot_df<- plot_df %>% mutate(scenario_num = case_when(
  scenario == "case1" ~ 1,
  scenario == "case2" ~ 2,
  scenario == "case3" ~ 3,
  scenario == "case4" ~ 4
))


p4<-plot_df %>% filter(parameter %in% c("pi1", "pi2")) %>%  ggplot(aes(x = factor(scenario_num), y = estimate, fill= factor(scenario_num))) +
  geom_boxplot() +
  geom_hline(aes(yintercept = true_value), 
       color = "red", linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ parameter, scales = "free_y",
      labeller = labeller(parameter = c(pi1 = "pi[1]", pi2 = "pi[2]"), .default = label_parsed)) +
  labs(x = "",y = "Estimate") +
  theme_bw() +
  theme(strip.text = element_text(size = 16),legend.position = "none")


p5<-plot_df %>% filter(parameter %in% c("lambda1", "lambda2")) %>%  ggplot(aes(x = factor(scenario_num), y = estimate, fill= factor(scenario_num))) +
  geom_boxplot() +
  geom_hline(aes(yintercept = true_value), 
             color = "red", linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ parameter, scales = "free_y",
             labeller = labeller(parameter = c(lambda1 = "lambda[1]", lambda2 = "lambda[2]"), .default = label_parsed)) +
  labs(x = "",y = "Estimate") +
  theme_bw() +
  theme(strip.text = element_text(size = 16),legend.position = "none")


#mu, component 1
plot_df_mu<-plot_df %>% filter(parameter %in% c("mu1_1","mu1_2", "mu1_3", "mu1_4", "mu1_5",
                                                "mu2_1","mu2_2", "mu2_3", "mu2_4","mu2_5")) 

#I want to re-organize by cluster
plot_df_mu <- plot_df_mu %>%
  mutate(Cluster = factor(paste("Cluster", substr(parameter, 3, 3))))

plot_df_mu <- plot_df_mu %>%
  mutate(
    # 'mu1_3' -> '3')
    mu_index = sub(".*_", "", parameter),
    # Create a plotmath-ready string (e.g., 'mu[3]')
    mu_label = paste0("mu[", mu_index, "]")
  )

p6<-plot_df_mu %>%  ggplot(aes(x = factor(scenario_num), y = estimate, fill= factor(scenario_num))) +
  geom_boxplot() +
  geom_hline(aes(yintercept = true_value), 
             color = "red", linetype = "dashed", size = 0.7) +
  facet_grid(Cluster ~ mu_label, scales = "free", 
             labeller = labeller(Cluster = label_value, mu_label = label_parsed)) +
  labs(x = "Scenario",y = "Estimate") +
  theme_bw() +
  theme( strip.text = element_text(size = 16),legend.position = "none")


jpeg(file = paste0(path, "Ombmm_pi_estimates_hp.jpeg"),  
    width = 800, 
    height = 400) 
p4
dev.off()

jpeg(file = paste0(path, "Ombmm_lambda_estimates_hp.jpeg"),  
     width = 800, 
     height = 400) 
p5
dev.off()


jpeg(file = paste0(path, "Ombmm_mu_estimates_hp.jpeg"),  
     width = 800, 
     height = 400) 
p6
dev.off()


jpeg(file = paste0(path, "Ombmm_estimates_patched_hp.jpeg"),  
     width = 800, 
     height = 800) 
((p4 + p5)/p6)
dev.off()



########################################################################################
########################################################################################

# Feature Saliency: PIP
results_list_pip <- lapply(files, function(f) {
  x <- readRDS(f)
  
  if (is.null(x$PIP)) return(NULL)
  
  list(
    scenario = x[[1]]$scenario,
    replicate = x[[1]]$replicate,
    pip = x$PIP
  )
})


pip_df <- bind_rows(lapply(results_list_pip, function(x) {
  if (is.null(x)) return(NULL)
  
  tibble(
    scenario = x$scenario,
    replicate = x$replicate,
    variable = paste0("V", seq_along(x$pip)),
    pip = x$pip
  )
}))


# Scenarios 1-4
pip_df <- pip_df %>% mutate(scenario_num = case_when(
  scenario == "case1" ~ 1,
  scenario == "case2" ~ 2,
  scenario == "case3" ~ 3,
  scenario == "case4" ~ 4
))


#Split by scen: make a function that outputs graphs for each scenario
# 
# pip_summary_plots<-function(df, scen){
#   
# #Long format of the data and extract pip--TO BE USED IN GGPLOT
#   df_scen <- df %>% filter(scenario_num == scen)
#   pip_scen<- as.data.frame(do.call(rbind, df_scen$pip))
#   pip_long_scen<- pivot_longer(pip_scen, cols = starts_with("V"),
#                                 names_to = "Var", values_to = "PIP")
#   
# 
# #- PLOTS-
#   #1. Lollipop
#   avgpip_scen<- pip_scen %>%
#     pivot_longer(everything(), names_to = "variable", values_to = "pip") %>%
#     summarise(
#       meanpip = mean(pip, na.rm = TRUE),
#       sdpip   = sd(pip, na.rm = TRUE),
#       .by = variable
#     )
#   
#   p1<-avgpip_scen %>% ggplot(aes(x=variable, y=meanpip)) +
#     geom_segment(aes(x=variable ,xend=variable, y=0, yend=meanpip), color="skyblue") +
#     geom_point(size=2,alpha=0.6, color="#69b3a2") +
#     coord_flip() +
#     theme_ipsum() +
#     theme(
#       panel.grid.minor.y = element_blank(),
#       panel.grid.major.y = element_blank(),
#       axis.text = element_text(size=10),
#       legend.position="none"
#     ) +
#     ylim(0,1) +
#     labs(
#       title = paste0("Scenario", " ", scen),
#       #subtitle = "Across simulation replicates",
#       x = "",
#       y = "")
#   
#   # 
#   # #2. Density (ridges)
#   # p2<- pip_long_scen %>%
#   #   ggplot(aes(x = PIP, y = factor(Var))) +
#   #   geom_density_ridges(
#   #     fill = "cadetblue3",
#   #     alpha = 0.7,
#   #     color = "cadetblue4"
#   #   ) +
#   #   labs(
#   #     #title = "Distribution of Posterior Inclusion Probabilities",
#   #     #subtitle = "Across simulation replicates",
#   #     x = "",
#   #     y = "Features") +
#   #   theme_ridges()
#   
#   return(p1)
# }
pip_summary_plots <- function(df, scen) {
  
  df_scen <- df %>% filter(scenario_num == scen)
  
  avgpip_scen <- df_scen %>%
    summarise(
      meanpip = mean(pip, na.rm = TRUE),
      sdpip   = sd(pip, na.rm = TRUE),
      .by = variable
    )
  avgpip_scen <- avgpip_scen %>%
    mutate(
      variable = factor(variable, levels = paste0("V", 1:10))
    )
  
  p1 <- avgpip_scen %>%
    ggplot(aes(x = variable, y = meanpip)) +
    geom_segment(
      aes(x = variable, xend = variable, y = 0, yend = meanpip),
      color = "skyblue"
    ) +
    geom_point(size = 2, alpha = 0.8, color = "#69b3a2") +
    coord_flip() +
    theme_ipsum() +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text = element_text(size = 10),
      legend.position = "none"
    ) +
    ylim(0, 1) +
    labs(
      title = paste0("Scenario ", scen),
      x = "",
      y = "Mean PIP"
    )
  
  return(p1)
}



case1<-pip_summary_plots(pip_df, scen = 1)
case2<-pip_summary_plots(pip_df, scen = 2)
case3<-pip_summary_plots(pip_df, scen = 3)
case4<-pip_summary_plots(pip_df, scen = 4)
  

#https://patchwork.data-imaginist.com/articles/guides/layout.html
p1a<- (case1 + case2)/ (case3 + case4) + plot_layout(axes = "collect") + plot_annotation(title = 'Distribution of Posterior Inclusion Probabilities',
                                                                                subtitle = "Across simulation replicates") + plot_layout(axes = "collect")


p2a<-case1 + case2 + plot_annotation(title = 'Distribution of Posterior Inclusion Probabilities',
                                               subtitle = "Across simulation replicates") + plot_layout(axes = "collect")

p2b<-case3+ case4 + plot_annotation(title = 'Distribution of Posterior Inclusion Probabilities',
                                               subtitle = "Across simulation replicates") + plot_layout(axes = "collect")


jpeg(file = paste0(path, "Ombmm_scen1_2_lol_hp.jpeg"),  
     width = 800, 
     height = 800) 
p1a
dev.off()



#Summary tables

pip_df1 <- pip_df %>%
  mutate(
    Var = lapply(pip, function(x) paste0("V", seq_along(x)))
  ) %>%
  unnest(c(Var, pip)) %>%
  rename(PIP = pip)



#False inclusion rate
pip_df %>%
  filter(
    scenario %in% c("case3", "case4"),
    variable %in% c("V5", "V6", "V7", "V8", "V9", "V10")
  ) %>%
  group_by(scenario) %>%
  summarise(
    FIR_05 = mean(pip > 0.5, na.rm = TRUE),
    FIR_09 = mean(pip > 0.9, na.rm = TRUE),
    max_pip = max(pip, na.rm = TRUE),
    mean_pip = mean(pip, na.rm = TRUE),
    .groups = "drop"
  )
# # A tibble: 1 × 4
# A tibble: 2 × 5
# scenario FIR_05 FIR_09 max_pip mean_pip
# <chr>     <dbl>  <dbl>   <dbl>    <dbl>
#   1 case3     0.170  0.168       1    0.177
# 2 case4     0.17   0.168       1    0.173


pip_df %>%
  filter(scenario %in% c("case3", "case4")) %>%
  group_by(variable) %>%
  summarise(mean_pip = mean(pip)) %>%
  arrange(desc(mean_pip))


#Facet by scenarios--just the density
# 
# pip_long <- est_info_pip %>%
#   mutate(
#     Var = lapply(pip, function(x) paste0("V", seq_along(x)))
#   ) %>%
#   unnest(c(Var, pip)) %>%
#   rename(PIP = pip)

# pip_long %>% ggplot(aes(x = PIP, y = Var)) +
#   geom_density_ridges(bandwidth = 0.1,
#     from = 0,
#     to = 1,
#     fill = "cadetblue3",
#     alpha = 0.7,
#     color = "cadetblue4") +
#   facet_wrap(~ scenario) +
#   scale_x_continuous(breaks = seq(0, 1, by = 0.25),expand = c(0, 0)) +
#   labs(
#     title = "Posterior Inclusion Probability Distributions",
#     subtitle = "Across simulation replicates",
#     x = "Posterior Inclusion Probability",
#     y = "Features") +
#   theme_ridges() +
#   theme(legend.position = "none",
#     strip.text = element_text(face = "bold")
#   )



######################TRACEPLOTS#############################################
source("R/helpers.R")
source("R/bayesmbmm_OFMM_SF_v2.R")
set.seed(2025)
datlist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_dpgs1_4.rds")
alloclist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_dpgs1_4.rds")
#contains tibble with true parameters

scenariosdat<- alloclist[[1]] #tibble of true parameters (mu, lambda, pi)
true_params1<- get_true_params(scenariosdat, "case2")
true_params2<- get_true_params(scenariosdat, "case4")

rid<-sample(1:100, 1)


start_time <- Sys.time()
#changing prior for lambda because we have different values
fit2<-Ombmm_FeatSaliency(
  x = datlist$case2[[rid]],
  K = 10,
  niter = 10000,
  nburnin = 3000,
  nchain = 1,
  init_list = NULL, 
  g =  5, h = 3, a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
  a_rho = 1, b_rho = 20, e0 = 0.01,
  priorOnAlpha = "no",
  p.true = true_params1$true_pi, plot_trace = TRUE)

fit4<-Ombmm_FeatSaliency(
  x = datlist$case4[[rid]],
  K = 10,
  niter = 10000,
  nburnin = 3000,
  nchain = 1,
  init_list = NULL, 
  g =  5, h = 3, a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
  a_rho = 1, b_rho = 20, e0 = 0.01,
  priorOnAlpha = "no",
  p.true = true_params2$true_pi, plot_trace = TRUE)

end_time <- Sys.time()
runtime_sec <- as.numeric(difftime(end_time, start_time, units = "secs")) # ~43 minutes


pdf(file = "simulations/traceplots/Case2_diagnostics.pdf")
#case 2
color_scheme_set("brightblue")
s_raw<-as.matrix(fit2$samples_keep)
s<-as.matrix(fit2$modeloutput_perchain$parameters.ECR.mcmc)

mcmc_trace(s_raw[,paste0("pi[", 1:10, "]")], 
           pars = paste0("pi[", 1:10, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for pi (all) case 2")

autocorr.plot(s_raw[,paste0("lambda[", 1:10, "]")])

mcmc_trace(s_raw[,paste0("lambda[", 1:10, "]")], 
           pars = paste0("lambda[", 1:10, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for lambda (all) case 2")

#Ordered and final
Kmode<-fit2$K
## plot trace of lambda
mcmc_trace(s[,paste0("lambda[", 1:Kmode, "]")], 
           pars = paste0("lambda[", 1:Kmode, "]"),
           facet_args = list(ncol = 1)) +
  ggtitle("Trace plots for lambda (2 clusters) case 2")


## plot trace of mu for feature 1

mu_samples1<-s[,paste0("mu[", 1:Kmode, ", ", 1, "]")]
color_scheme_set("brightblue")
mcmc_trace(mu_samples1, 
           pars = paste0("mu[", 1:Kmode, ", ", 1, "]"),
           facet_args = list(ncol = 1)) +
  ggtitle("Trace plots for mu, feature 1, case 2")


## plot trace of mu for feature 2

mu_samples2<-s[,paste0("mu[", 1:Kmode, ", ", 2, "]")]
color_scheme_set("brightblue")
mcmc_trace(mu_samples2, 
           pars = paste0("mu[", 1:Kmode, ", ", 2, "]"),
           facet_args = list(ncol = 1)) +
  ggtitle("Trace plots for mu, feature 2, case 2")

dev.off()













pdf(file = "simulations/traceplots/Case4_diagnostics.pdf")
#case 4
color_scheme_set("brightblue")
s_raw<-as.matrix(fit4$samples_keep)
s<-as.matrix(fit4$modeloutput_perchain$parameters.ECR.mcmc)


mcmc_trace(s_raw[,paste0("pi[", 1:10, "]")], 
           pars = paste0("pi[", 1:10, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for pi (all) case 4")

autocorr.plot(s_raw[,paste0("lambda[", 1:10, "]")])

#Iterations with non-empty components

mcmc_trace(s_raw[,paste0("lambda[", 1:10, "]")], 
           pars = paste0("lambda[", 1:10, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for lambda (all) case 4")

#Ordered and final
Kmode<-fit4$K
## plot trace of lambda
mcmc_trace(s[,paste0("lambda[", 1:Kmode, "]")], 
           pars = paste0("lambda[", 1:Kmode, "]"),
           facet_args = list(ncol = 1)) +
  ggtitle("Trace plots for lambda (2 clusters) case 4")



## plot trace of mu for feature 1

mu_samples1<-s[,paste0("mu[", 1:Kmode, ", ", 1, "]")]
color_scheme_set("brightblue")
mcmc_trace(mu_samples1, 
           pars = paste0("mu[", 1:Kmode, ", ", 1, "]"),
           facet_args = list(ncol = 1)) +
  ggtitle("Trace plots for mu, feature 1, case 4")


## plot trace of mu for feature 2

mu_samples2<-s[,paste0("mu[", 1:Kmode, ", ", 2, "]")]
color_scheme_set("brightblue")
mcmc_trace(mu_samples2, 
           pars = paste0("mu[", 1:Kmode, ", ", 2, "]"),
           facet_args = list(ncol = 1)) +
  ggtitle("Trace plots for mu, feature 2, case 4")

dev.off()




