#In this script, we combine all outputs after cluster run and summarize results
# Programmer: Carmen Rodriguez C.
# Date updated: 3/17/2026
#=========================
library(nimble)
library(ggplot2)
library(GGally)
library(tidyverse)
library(purrr)
library(patchwork)
library(coda)
library(bayesplot)
library(parallel)
library(hrbrthemes)


#NOTE: lambda_est and mu_est, need to check label switching... i used pi_true as anchor but they're all equal

#All reps
#------First we look at model parameters for all
files_path<-"/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs"
files <- list.files(files_path, full.names = TRUE, pattern = "^SAMBOACSmockFinal_.*\\.rds$")
test<-readRDS(files[1])

#check file size
#i saved the samples for further exploring 
sizes <- file.info(files)$size
summary(sizes)

###########################################
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

res_list <- lapply(files, safe_read)
res_df <- do.call(rbind, res_list)

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
      PIP = x$PIP
    )
    
  }, error = function(e) {
    message("FAILED: ", f)
    return(NULL)
  })
}

pip_list <- lapply(files, extract_pip_file)

# Remove NULLs
pip_list <- Filter(Negate(is.null), pip_list)

# Combine
pip_df <- do.call(rbind, pip_list)

#saveRDS(res_df, file.path("~/Mbeta_Project/Sim_results", "SAMBOacsmock_res_df_03.26.rds"))
#saveRDS(pip_df, file.path("~/Mbeta_Project/Sim_results", "SAMBOacsmock_pip_df_03.26.rds"))



#-----Read in files----
#Going through the results to see Kmode, ARI,
res_df<-readRDS("~/Mbeta_Project/Sim_results/SAMBOacsmock_res_df_03.26.rds")

#quick check
res_df %>%
  group_by(scen) %>%
  summarise(
    n = n(),
    failures = sum(failed, na.rm = TRUE),
    ARI_na = sum(is.na(ARI)),
    MSE_mu_na = sum(is.na(MSE_mu))
  )

#10 reps failed

#Check ARI
ggplot(res_df %>% filter(!failed), aes(x = scen, y = ARI)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Clustering Performance (ARI)")

#MSE for the mean
#not bad across-- coud be better
ggplot(res_df %>% filter(!failed), aes(x = scen, y = MSE_mu)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Parameter Recovery (MSE_mu)")


#Runtime
ggplot(res_df %>% filter(!failed), aes(x = scen, y = (runtime/3600))) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Runtime by Scenario (hours)")


summary_table <- res_df %>%
  filter(!failed) %>%
  group_by(scen) %>%
  summarise(
    ARI_mean = mean(ARI, na.rm = TRUE),
    ARI_sd   = sd(ARI, na.rm = TRUE),
    
    MSE_mu_mean = mean(MSE_mu, na.rm = TRUE),
    MSE_mu_sd   = sd(MSE_mu, na.rm = TRUE),
    
    K_mean = mean(Kmode, na.rm = TRUE),
    
    runtime_mean = mean(runtime, na.rm = TRUE)
  )

res_df <- res_df %>%
  left_join(
    scenariosdat %>%
      select(dgp_id, p, overlap),
    by = c("scen" = "dgp_id")
  )





#PIP exploration
#Ground truth
pip_truth_list <- list(
  MockACS1 = rep(1, 12),
  MockACS2 = rep(1, 12),
  MockACS3 = rep(1, 12),
  MockACS4 = c(rep(1, 12), rep(0, 8)),
  MockACS5 = c(rep(1, 12), rep(0, 8))
)

pip_df$pip_truth <- mapply(
  function(s, f) pip_truth_list[[s]][f],
  pip_df$scen,
  pip_df$feature
)


#Scenarios 1-3 (all relevant)
pip_df %>%
  filter(scen %in% c("MockACS1","MockACS2","MockACS3")) %>%
  group_by(scen) %>%
  summarise(mean_PIP = mean(PIP),
            min_PIP = min(PIP),
            sd_PIP = sd(PIP))


#Scenarios 4-5 (8 noise variables)

ggplot(pip_df %>% filter(scen %in% c("MockACS4","MockACS5")),
  aes(x = PIP, fill = factor(pip_truth))) +
  geom_density(alpha = 0.4) +
  facet_wrap(~scen) +
  labs(fill = "Relevant")

pip_df %>%
  filter(pip_truth == 0) %>%
  group_by(scen) %>%
  summarise(mean_noise_PIP = mean(PIP))

# # A tibble: 2 × 2
# scen     mean_noise_PIP
# <chr>             <dbl>
#   1 MockACS4       0.00774 
# 2 MockACS5       0.000729

pip_df %>%
  filter(pip_truth == 1) %>%
  group_by(scen) %>%
  summarise(mean_PIP = mean(PIP))
# scen     mean_PIP
# <chr>       <dbl>
#   1 MockACS1    1.00 
# 2 MockACS2    1.00 
# 3 MockACS3    1.00 
# 4 MockACS4    0.997
# 5 MockACS5    1    



pip_df %>%
  filter(scen %in% c("MockACS4","MockACS5")) %>%
  group_by(scen, pip_truth) %>%
  summarise(
    mean_PIP = mean(PIP),
    sd_PIP = sd(PIP)
  )
# scen     pip_truth mean_PIP sd_PIP
# <chr>        <dbl>    <dbl>  <dbl>
#   1 MockACS4         0 0.00774  0.0648
# 2 MockACS4         1 0.997    0.0333
# 3 MockACS5         0 0.000729 0.0179
# 4 MockACS5         1 1        0     



avgpip_scen4 <- pip_df %>%
  filter(scen == "MockACS4") %>%
  summarise(
    meanpip = mean(PIP, na.rm = TRUE),
    sdpip   = sd(PIP, na.rm = TRUE),
    .by = feature
  )

avgpip_scen5 <- pip_df %>%
  filter(scen == "MockACS5") %>%
  summarise(
    meanpip = mean(PIP, na.rm = TRUE),
    sdpip   = sd(PIP, na.rm = TRUE),
    .by = feature
  )
p1 <- avgpip_scen4 %>% ggplot(aes(x=feature, y=meanpip)) +
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
    title = paste0("Scenario", " ","MockACS4" ),
    #subtitle = "Across simulation replicates",
    x = "",
    y = "")



p2 <- avgpip_scen5 %>% ggplot(aes(x=feature, y=meanpip)) +
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
    title = paste0("Scenario", " ", "MockACS5"),
    #subtitle = "Across simulation replicates",
    x = "",
    y = "")


plotspip<- (p1 +p2) + plot_layout(axes = "collect") + plot_annotation(title = 'Distribution of Posterior Inclusion Probabilities',
                                                                                                             subtitle = "Across simulation replicates") + plot_layout(axes = "collect")
path<-"~/bayesmbmm/simulations/figures/"
jpeg(file = paste0(path, "sambo_mockacs_scen4_5_lol.jpeg"),  
     width = 800, 
     height = 800) 
plotspip
dev.off()



#Look at a few trace plots, especially for last 2 data cases
results_subset<- lapply(files[c(327, 400,450)], function(f) readRDS(f))
s1<-results_subset[[3]]$samplesall
K<-20

## plot trace of pi
pi_samples<-s1[,paste0("pi[", 1:K, "]")]
color_scheme_set("mix-brightblue-gray")
mcmc_trace(pi_samples, 
           pars = paste0("pi[", 1:10, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for pi")

mcmc_trace(pi_samples, 
           pars = paste0("pi[", 11:20, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for pi")
#2,4,10, 15, 17

effectiveSize(pi_samples)
## plot trace of mu for feature

mu_samples1<-s1[,paste0("mu[", 2, ", ", 1:12, "]")]
color_scheme_set("brightblue")
mcmc_trace(mu_samples1, 
           pars = paste0("mu[", 2, ", ", 1:12, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for mu, comp 2")

effectiveSize(mu_samples1) #good!


# all around 0.5
mu_samples2<-s1[,paste0("mu[", 2, ", ", 13:20, "]")]
mcmc_trace(mu_samples2, 
           pars = paste0("mu[", 2, ", ", 13:20, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for mu, comp 2")



#lambda
lambda_samples<-s1[,paste0("lambda[", 1:K, "]")]

mcmc_trace(lambda_samples, 
           pars = paste0("lambda[", 1:10, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for lambda")

#2,4,10

mcmc_trace(lambda_samples, 
           pars = paste0("lambda[", 11:20, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for lambda")


#15, 17

samplesSummary(lambda_samples)


#Look at the estimated values vs. truth
# > scenariosdat[5,"lambda_vec"][[1]]
# [[1]]
# [1] 10 18 25 35 60
# > samplesSummary(lambda_samples)
# Mean   Median   St.Dev. 95%CI_low 95%CI_upp
# lambda[2]  10.09730 10.09417 0.1656579  9.781517  10.43216
# lambda[4]  24.43130 24.43401 0.3929528 23.657829  25.20441
# lambda[10] 57.86652 57.85442 1.0218273 55.874893  59.93991
# lambda[15] 17.67509 17.66989 0.2976107 17.110082  18.27935
# lambda[17] 34.34708 34.34162 0.5455622 33.272030  35.42016
# 

lambda_est<-results_subset[[3]]$lambda_est
#2 and 4 are swapped

#y<-c(10.09730, 17.67509, 24.43130, 34.34708, 57.86652)




# #Failed jobs
# failed_files <- files[sizes <= 1e6]  # from earlier
# 
# parse_info <- function(f) {
#   fname <- basename(f)
#   
#   # Extract scenario and rep using regex
#   scen <- sub(".*SAMBOACSmockFinal_(MockACS[0-9]+)_rep[0-9]+\\.rds", "\\1", fname)
#   rep  <- as.integer(sub(".*_rep([0-9]+)\\.rds", "\\1", fname))
#   
#   data.frame(scen = scen, rep = rep)
# }
# 
# failed_jobs <- do.call(rbind, lapply(failed_files, parse_info))




#TEST OF ONLY 3 REPS
#------First we look at model parameters for all
files_path<-"/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs"

datlist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/MockACSUpdated/datlistAll_mockacsdpgs1_5.rds") 
alloclist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/MockACSUpdated/alloclistAll_mockacsdpgs1_5.rds")
# Define all scenario × replicate pairs
scenariosdat<- alloclist[[1]] #tibble of true parameters (mu, lambda, pi)

scenarios <- names(datlist)  
jobs <- expand.grid(scen = scenarios, rep = 1:3, stringsAsFactors = FALSE)


results <- lapply(seq_len(nrow(jobs)), function(i) {
  cid <- jobs$scen[i]
  rid <- jobs$rep[i]
  
  file_path <- sprintf(
    "%s/SAMBOACSmock_%s_rep%d.rds",
    files_path, cid, rid
  )
  
  if (file.exists(file_path)) {
    res <- readRDS(file_path)
    
    # attach metadata (IMPORTANT)
    res$scenario_name <- cid
    res$replicate <- rid
    
    return(res)
  } else {
    warning("Missing file: ", file_path)
    return(NULL)
  }
})


test<-results[[1]]
s<-as.matrix(test$samplesall)
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
#3, 6,8,11,16,

effectiveSize(pi_samples)
## plot trace of mu for feature

mu_samples1<-s[,paste0("mu[", 3, ", ", 1:12, "]")]
color_scheme_set("brightblue")
mcmc_trace(mu_samples1, 
           pars = paste0("mu[", 3, ", ", 1:12, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for mu, comp 3")

effectiveSize(mu_samples1) #very low!


#perfect-- all around 0.5
mu_samples2<-s[,paste0("mu[", 1, ", ", 1:12, "]")]
color_scheme_set("brightblue")
mcmc_trace(mu_samples2, 
           pars = paste0("mu[", 1, ", ", 1:12, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for mu, comp 1")



test2<-results[[4]] #p=20
s<-as.matrix(test2$samplesall)
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
#2,9,11,15,17


## plot trace of mu for feature -- burn-in to 5000

mu_samples1<-s[,paste0("mu[", 2, ", ", 1:12, "]")]
color_scheme_set("brightblue")
mcmc_trace(mu_samples1, 
           pars = paste0("mu[", 2, ", ", 1:12, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for mu, comp 2")


#identical- nonimportant vars
#Non-important variables (~0.5) behave nicely
mu_samples2<-s[,paste0("mu[", 2, ", ", 13:20, "]")]
color_scheme_set("brightblue")
mcmc_trace(mu_samples2, 
           pars = paste0("mu[", 2, ", ", 13:20, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for mu, comp 2, nonimportant vars")


#perfect-- all around 0.5
mu_samples2<-s[,paste0("mu[", 1, ", ", 1:12, "]")]
color_scheme_set("brightblue")
mcmc_trace(mu_samples2, 
           pars = paste0("mu[", 1, ", ", 1:12, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for mu, comp 1")




 
# summarize_results <- function(res_list) {
#   do.call(rbind, lapply(res_list, function(x) {
#     
#     data.frame(
#       scenario = x$scenario_name,
#       rep = x$replicate,
#       status = x$status,
#       Kmode = ifelse(x$status == "success", x$Kmode, NA),
#       runtime_sec = x$runtime_sec,
#       n_clusters = ifelse(x$status == "success", length(x$pi_est), NA),
#       mean_lambda = ifelse(x$status == "success", mean(x$lambda_est), NA),
#       mean_PIP = ifelse(x$status == "success", mean(x$PIP), NA),
#       stringsAsFactors = FALSE
#     )
#     
#   }))
# }
# 
# df_summary <- summarize_results(results)
#

###Other checks
df_runtime <- data.frame(
  scenario = sapply(results, function(x) x$scenario_name),
  rep = sapply(results, function(x) x$replicate),
  runtime = sapply(results, function(x) x$runtime))
  
pi_ess<-effectiveSize(results[[1]]$samplesall[, paste0("pi[", 1:K, "]")])
#Active clusters: 3, 6,8,11,16,

effectiveSize(results[[1]]$samplesall[, paste0("mu[", 3, ", ", 1:12, "]")]) 
effectiveSize(results[[1]]$samplesall[, paste0("mu[", 6, ", ", 1:12, "]")]) 
effectiveSize(results[[1]]$samplesall[, paste0("mu[", 8, ", ", 1:12, "]")]) 
effectiveSize(results[[1]]$samplesall[, paste0("mu[", 11, ", ", 1:12, "]")]) 
effectiveSize(results[[1]]$samplesall[, paste0("mu[", 16, ", ", 1:12, "]")]) 

check6<-results[[1]]$samplesall[, paste0("mu[", 6, ", ", 1:12, "]")]
samplesSummary(check6)






