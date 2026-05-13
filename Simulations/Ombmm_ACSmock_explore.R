#In this script, we combine all outputs after cluster run and summarize results
# Programmer: Carmen Rodriguez C.
# Date updated: 2/7/2026
#=========================
library(nimble)
library(ggplot2)
library(tidyverse)
library(bayesplot)
library(hrbrthemes)
#------First we look at model parameters for all
files_path<-"/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs"
files1 <- list.files(files_path, full.names = TRUE, pattern = "^Ombmmcases_old_case5_rep[0-9]+\\.rds$")
files2<-list.files(files_path, full.names = TRUE, pattern = "^Ombmmcases_old_case7_rep[0-9]+\\.rds$")

#test<-readRDS("/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/Ombmmcases_old_case5_rep1.rds")


safe_read <- function(f) {
  tryCatch({
    x <- readRDS(f)
    
    data.frame(
      scen = x$cid,
      rep = x$rid,
      runtime = x$runtime,
      Kmode = x$Kmode,
      mse_lambda = x$mse_lambda,
      MSE_mu = x$MSE_mu,
      ARI = x$ARI,
      failed = isTRUE(x$failed),
      file = basename(f)
    )
    
  }, error = function(e) {
    message("FAILED: ", f)
    
    return(data.frame(
      scen = NA,
      rep = NA,
      runtime = NA,
      Kmode = NA,
      mse_lambda = NA,
      MSE_mu = NA,
      ARI = NA,
      failed = TRUE,
      file = basename(f)
    ))
  })
}

res5_list <- lapply(files1, safe_read)
res7_list <- lapply(files2, safe_read)
res5_df <- do.call(rbind, res5_list)
res7_df <- do.call(rbind, res7_list)

res_all<-rbind(res5_df,res7_df)
res_all$K_true<-5

res_all %>%
  group_by(scen) %>%
  summarise(
    mean_ARI = mean(ARI),
    sd_ARI   = sd(ARI),
    mean_K   = mean(Kmode),
    P_correct_K = mean(Kmode == K_true),
    P_over  = mean(Kmode > K_true),
    P_under = mean(Kmode < K_true))


ggplot(res_all, aes(x = ARI, fill = scen)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  facet_wrap(~ scen) +
  theme_minimal()


ggplot(res_all, aes(x = Kmode)) +
  geom_bar() +
  facet_wrap(~ scen) +
  geom_vline(
    aes(xintercept = K_true),
    linetype = "dashed",
    color = "red"
  ) +
  theme_minimal()

#MSE for the mean
#not bad across-- coud be better
ggplot(res_all %>% filter(!failed), aes(x = scen, y = MSE_mu)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Parameter Recovery (MSE_mu)")

#Runtime
ggplot(res_all %>% filter(!failed), aes(x = scen, y = (runtime/3600))) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Runtime by Scenario (hours)")




#PIP
##############################################
extract_pip_file <- function(f) {
  tryCatch({
    x <- readRDS(f)
    # Skip failed or missing PIP
    if (is.null(x$PIP)) return(NULL)
    
    data.frame(
      scen = x$cid,
      rep = x$rid,
      feature = seq_along(x$PIP),
      PIP = x$PIP,
      Kmode = x$Kmode,
      mse_lambda = x$mse_lambda,
      MSE_mu = x$MSE_mu
    )
    
  }, error = function(e) {
    message("FAILED: ", f)
    return(NULL)
  })
}

pip_list5 <- lapply(files1, extract_pip_file)
pip_list7 <- lapply(files2, extract_pip_file)


# Combine
pip_df5 <- do.call(rbind, pip_list5)
pip_df7 <- do.call(rbind, pip_list7)
pip_df<-rbind(pip_df5,pip_df7)


#PIP exploration
# #Ground truth--wrong
# pip_truth_list <- list(
#   case5 = rep(1, 12),
#   case7 = c(rep(1, 12), rep(0, 8))
# )
# 
# pip_df$pip_truth <- mapply(
#   function(s, f) pip_truth_list[[s]][f],
#   pip_df$scen,
#   pip_df$feature
# )

#mean pip for each feature
pip_df %>%
  filter(scen == "case5") %>%
  group_by(feature) %>% 
  summarise(mean_PIP = mean(PIP),
            min_PIP = min(PIP),
            sd_PIP = sd(PIP))



pip_df %>%
  filter(scen == "case7") %>%
  group_by(feature) %>% 
  summarise(mean_PIP = mean(PIP),
            min_PIP = min(PIP),
            sd_PIP = sd(PIP))



avgpip_scen7 <- pip_df %>%
  filter(scen == "case7") %>%
  summarise(
    meanpip = mean(PIP, na.rm = TRUE),
    sdpip   = sd(PIP, na.rm = TRUE),
    .by = feature
  )

plot1<- avgpip_scen7 %>% ggplot(aes(x=feature, y=meanpip)) +
  geom_segment(aes(x=feature ,xend=feature, y=0, yend=meanpip), color="skyblue") +
  geom_point(size=2,alpha=0.6, color="#69b3a2") +
  coord_flip() +
  theme_ipsum() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(size=10),
    legend.position="none"
  ) +
  ylim(0,1) +
  labs(
    title = paste0("Scenario", " ", 6),
    #subtitle = "Across simulation replicates",
    x = "",
    y = "")


png(filename = "~/bayesmbmm/simulations/figures/lollipop_scenario6.png")
plot1
dev.off()




avgpip_scen5 <- pip_df %>%
  filter(scen == "case5") %>%
  summarise(
    meanpip = mean(PIP, na.rm = TRUE),
    sdpip   = sd(PIP, na.rm = TRUE),
    .by = feature
  )

avgpip_scen5 %>% ggplot(aes(x= feature, y=meanpip)) +
  geom_segment(aes(x=feature ,xend=feature, y=0, yend=meanpip), color="skyblue") +
  geom_point(size=2,alpha=0.6, color="#69b3a2") +
  coord_flip() +
  theme_ipsum() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(size=10),
    legend.position="none"
  ) +
  ylim(0,1) +
  labs(
    title = paste0("Scenario", " ", 5),
    #subtitle = "Across simulation replicates",
    x = "",
    y = "")


######Filter when Kmode = Ktrue
pip_df_5<- pip_df %>% filter(Kmode == 5)

pip_df_5 %>%
  filter(scen == "case5") %>%
  group_by(feature) %>% 
  summarise(mean_PIP = mean(PIP),
            min_PIP = min(PIP),
            sd_PIP = sd(PIP))
pip_df_5 %>%
  filter(scen == "case7") %>%
  group_by(feature) %>% 
  summarise(mean_PIP = mean(PIP),
            min_PIP = min(PIP),
            sd_PIP = sd(PIP))


#MSE for the mean
res_all_5<- res_all %>% filter(Kmode == 5)
reps_5<- res_all_5$rep

scenario5_6_K5<-ggplot(res_all_5 %>% filter(!failed), aes(x = scen, y = MSE_mu)) +
  geom_boxplot() +
  facet_wrap(~ scen) +
  theme_minimal() +
  labs(title = "Parameter Recovery (MSE mu)")



ggplot(res_all_5 %>% filter(!failed), aes(x = scen, y = mse_lambda)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Parameter Recovery (MSE lambda)")


ggplot(res_all_5, aes(x = ARI, fill = scen)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  facet_wrap(~ scen) +
  theme_minimal()






#SAVE RESULTS
saveRDS(list(res_all, pip_df), file = "/n/home03/crodriguezcabrera/Mbeta_Project/Sim_results/SAMBO_acsmock_case5_7_3.20.26.rds")


datlist5<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistAll_case5.rds")
alloclist5<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case5.rds")


#5
dat5_11<-as.data.frame(datlist5[[11]])
#4
dat5_15<-as.data.frame(datlist5[[15]])
#3
dat5_23<-as.data.frame(datlist5[[23]])

alloc_11<-factor(alloclist5[[11]])
alloc_15<-factor(alloclist5[[15]])
alloc_23<-factor(alloclist5[[23]])


p1<-ggpairs(dat5_11,
            ggplot2::aes(color = alloc_11, alpha = 0.5),
            upper = list(continuous = "blank"),
            lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 10)

p2<-ggpairs(dat5_15,
            ggplot2::aes(color = alloc_15, alpha = 0.5),
            upper = list(continuous = "blank"),
            lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 10)

p3<-ggpairs(dat5_23,
                ggplot2::aes(color = alloc_23, alpha = 0.5),
              upper = list(continuous = "blank"),
                lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 10)

pdf(file ="data/5clustACS_mock_same_4.22.26.pdf",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 10) # The height of the plot in inches
p1
p2
p3
dev.off()








## Traceplots
check<-readRDS(files1[1])
s<-check$samplesall
K<-20

## plot trace of pi
pi_samples<-s[,paste0("pi[", 1:K, "]")]
color_scheme_set("mix-brightblue-gray")
mcmc_trace(pi_samples, 
           pars = paste0("pi[", 1:10, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for pi")

mcmc_trace(pi_samples, 
           pars = paste0("pi[", 11:20, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for pi")



lambda_samples<-s[,paste0("lambda[", 1:K, "]")]
color_scheme_set("brightblue")
mcmc_trace(lambda_samples, 
           pars = paste0("lambda[", 1:10, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for lambda")
mcmc_trace(lambda_samples, 
           pars = paste0("lambda[", 11:20, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for lambda")

#choosing comp 12 and 18-- worst trace plots for lambda--why?
#estimates for lambda are really bad!

samplesSummary(lambda_samples)

## plot trace of mu for feature

mu_samples1<-s[,paste0("mu[", 12, ", ", 1:12, "]")]
color_scheme_set("brightblue")
mcmc_trace(mu_samples1, 
           pars = paste0("mu[", 12, ", ", 1:12, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for mu, comp 12")



mu_samples2<-s[,paste0("mu[", 18, ", ", 1:12, "]")]
color_scheme_set("brightblue")
mcmc_trace(mu_samples2, 
           pars = paste0("mu[", 18, ", ", 1:12, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for mu, comp 18")





### Test with new DGP--2/24/26 ### See ~/FASRC/testing.R
acs1<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistmockACS1.rds") 
acs2<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/datlistmockACS2.rds")

#only 6/8 jobs completed
files1 <- list.files(files_path, full.names = TRUE, pattern = "^OmbmmACSmock_testing_.*\\.rds$")

test2<-readRDS("/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/OmbmmACSmock_testing_job4.rds")
str(test2)
###JOB4: concentration parameter 
# n<-2000
# Ktrue <-5
# p<-12 #important vars
# nk<-c(600, 500, 400, 300, 200)
# pi_true1<- nk/n #0.30 0.25 0.20 0.15 0.10
# lambda_vec2<- rep(25,5)

test3<-readRDS("/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/OmbmmACSmock_testing_job8.rds")

ggpairs(acs1[[4]], columns = 1:12,
        ggplot2::aes(alpha = 0.5),
        upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 10)


set.seed(2025 + 4)
acs1_s4<-Ombmm_FeatSaliency(
    x = acs1[[4]],
    K = 20,
    niter = 15000,
    nburnin = 5000,
    nchain = 1,
    init_list = NULL,
    g = 5, h = 3,
    a_mu0 = 1, b_mu0 = 1,
    a_mu1 = 2, b_mu1 = 2,
    a_rho = 1, b_rho = 20,
    e0 = 0.01,p.true = c(0.30, 0.25, 0.20, 0.15, 0.10),
    priorOnAlpha = "no",
    plot_trace = TRUE
  )


#x<- apply(acs1[[4]], 2, qlogis)
#fitclust<-Mclust(x, G = 1:20)


acs2_s4<-Ombmm_FeatSaliency(
  x = acs2[[4]],
  K = 20,
  niter = 15000,
  nburnin = 5000,
  nchain = 1,
  init_list = NULL,
  g = 5, h = 3,
  a_mu0 = 1, b_mu0 = 1,
  a_mu1 = 2, b_mu1 = 2,
  a_rho = 1, b_rho = 20,
  e0 = 0.01,p.true = c(0.30, 0.25, 0.20, 0.15, 0.10),
  priorOnAlpha = "no",
  plot_trace = TRUE
)






#------First we look at model parameters for all
files_path<-"/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs"
files <- list.files(files_path, full.names = TRUE, pattern = "^OmbmmresultACSmock_cniter2_.*\\.rds$")


results_list <- lapply(files, function(f) readRDS(f))








