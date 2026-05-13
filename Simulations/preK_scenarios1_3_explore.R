
#sacct -j 47868422
#=================================================================
# Used simulated data from Cases 1-3 (see excel file with values)
#For each case,  we generated 100 datasets and calculated MSE, MAB and ARI for all model parameters
#In this script, we combine all outputs after cluster run and summarize results
# Programmer: Carmen Rodriguez C.
# Date updated: 2025/11/24
#===========================================================
library(ggplot2)
library(tidyverse)

files_path<-"/n/netscratch/stephenson_lab/Lab/crodriguez/MBeta_outputs"
files <- list.files(files_path, full.names = TRUE, pattern = "rds$")
results_list <- lapply(files, function(f) readRDS(f)[[1]])#Only keep the data.frames
results_all <- do.call(rbind, results_list)

str(results_all)

summary_df <- results_all %>%
  group_by(scenario) %>%
  summarise(across(c(mse_pi, mse_mu, mse_lambda, mab_pi, mab_mu, mab_lambda, ARIval),
                   list(mean = mean, sd = sd), .names = "{.col}_{.fn}"))


mse_long <- pivot_longer(results_all, cols = starts_with("mse_"),
                                names_to = "parameter", values_to = "MSE")

ggplot(mse_long, aes(x = parameter, y = MSE, fill = scenario)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA, position = position_dodge(0.9)) +
  scale_y_log10() +
  labs(title = "MSE across parameters and scenarios", y = "MSE (log scale)")



mab_long <- pivot_longer(results_all, cols = starts_with("mab_"),
                         names_to = "parameter", values_to = "MAB")

ggplot(mab_long, aes(x = parameter, y = MAB, fill = scenario)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA, position = position_dodge(0.9)) +
  scale_y_log10() +
  labs(title = "Absolute bias across parameters and scenarios", y = "MAB (log scale)")




#TO BE CHECKED

ggplot(results_all[results_all$scenario == "case1", ],
       aes(x = factor(replica), y = pi_est)) +
  geom_boxplot() +
  geom_hline(yintercept = pi_true, linetype = "dashed") +
  labs(title = "Parameter recovery per replicate (Scenario 1)",
       x = "Replicate", y = "Estimated π")

