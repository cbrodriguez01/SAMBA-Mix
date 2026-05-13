#SAMBO stress test
#In this script, we combine all outputs after cluster run and summarize results
# Programmer: Carmen Rodriguez C.
# Date updated: 3/21/2026
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

#All reps
#------First we look at model parameters for all
files_path<-"/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs"
files <- list.files(files_path, full.names = TRUE, pattern = "^SAMBO_stress2_.*\\.rds$")


fits <- lapply(files, readRDS)

n_runs <- length(fits)

summary_df <- data.frame(
  rep = 1:n_runs,
  K = sapply(fits, function(f) f$Kmode),
  ARI = sapply(fits, function(f) f$ARI),
  runtime = sapply(fits, function(f) f$runtime)
)


hist(summary_df$K, main = "Estimated K")
hist(summary_df$ARI, main = "ARI")


summary_stats <- data.frame(
  K_mean = mean(summary_df$K),
  K_sd = sd(summary_df$K),
  ARI_mean = mean(summary_df$ARI),
  ARI_sd = sd(summary_df$ARI),
  runtime_mean = mean(summary_df$runtime)
)

summary_stats

pip_mat <- do.call(rbind, lapply(fits, function(f) f$PIP))
pip_mean <- colMeans(pip_mat)


pip_groups <- list(
  strong = 1:6,
  weak = 7:12,
  corr_noise = 13:16,
  pure_noise = 17:20
)

pip_group_summary <- sapply(pip_groups, function(idx) {
  mean(pip_mean[idx])
})

pip_group_summary



