library(tidyverse)
library(ggplot2)

#alloclist5<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case5.rds")
#alloclist7<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case7.rds")


alloclist<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/MockACSUpdated/alloclistAll_mockacsdpgs1_5.rds")
true_allocation<-alloclist[[2]]


dir <- "/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/"
files <- list.files(dir, pattern = "clustvarsel_mockacs_rep", full.names = TRUE)
files <- files[order(files)]

fits <- lapply(files, readRDS)

clustvarsel_fits3 <- lapply(fits, `[[`, "mockacs3")
clustvarsel_fits5 <- lapply(fits, `[[`, "mockacs5")

true_vars_case3 <- paste0("X", 1:12)
true_vars_case5 <- paste0("X", 1:12)
noise_vars_case5 <- paste0("X", 13:20)


str(clustvarsel_fits3[[1]])


summarize_output<-function(fits, true_alloc,sid){
  map_dfr(seq_along(fits), function(r){
    fit <- fits[[r]]
    z_true <- true_alloc[[r]]
    z_hat  <- fit$classification
    
    tibble(case = sid,
           rep = r,
           ARI = adjustedRandIndex(z_true, z_hat),
           K_hat = fit$G_hat,
           model= fit$modelName,
           subsetLength = length(fit$subset)) #to see if all variables were selected from the model
  })
}


res_case3 <- summarize_output(
  fits = clustvarsel_fits3,
  true_alloc = true_allocation$MockACS3,
  sid = "MockACS3"
)

res_case5 <- summarize_output(
  fits = clustvarsel_fits5,
  true_alloc = true_allocation$MockACS5,
  sid ="MockACS5"
)

res_all <- bind_rows(res_case3, res_case5)
res_all$K_true<-5
res_all$method<- "clustvarsel"


res_all %>%
  group_by(case) %>%
  summarise(
    mean_ARI = mean(ARI),
    sd_ARI   = sd(ARI),
    mean_K   = mean(K_hat),
    P_correct_K = mean(K_hat == K_true),
    P_over  = mean(K_hat > K_true),
    P_under = mean(K_hat < K_true),
    SubsetNum = mean(subsetLength))


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


ggplot(res_all, aes(x = subsetLength)) +
  geom_bar() +
  facet_wrap(~ case) +
  geom_vline(
    aes(xintercept = 12),
    linetype = "dashed",
    color = "red"
  ) +
  theme_minimal()



var_freq_case5 <- map_dfr(clustvarsel_fits5, function(fit){
  tibble(var = names(fit$subset))
}) %>%
  count(var) %>%
  mutate(freq = n / length(clustvarsel_fits5),
         type = ifelse(var %in% noise_vars_case5, "noise", "signal"))

ggplot(var_freq_case5,
       aes(x = reorder(var, freq), y = freq, fill = type)) +
  geom_col() +
  coord_flip() +
  theme_minimal()


#Save results to home folder
saveRDS(res_all,file = "/n/home03/crodriguezcabrera/Mbeta_Project/Sim_results/clustvarsel_sim_new.rds")





