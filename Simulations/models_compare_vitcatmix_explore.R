# VICatMix variable selection for tertiles of the data
library(VICatMix)

alloclist5<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case5.rds")
alloclist7<-readRDS("/n/home03/crodriguezcabrera/Mbeta_Project/Data_for_sims/alloclistAll_case7.rds")

true_allocation<-alloclist[[2]]

files_path<-"/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/"
files1 <- list.files(files_path,full.names = TRUE,
  pattern = "^VICatMix_cases5_7_tertilescase5_rep[0-9]+\\.rds$")

files2<- list.files(files_path,full.names = TRUE,
  pattern = "^VICatMix_cases5_7_tertilescase7_rep[0-9]+\\.rds$"
)


results_list_case5 <- lapply(files1, function(f) readRDS(f)) 
results_list_case7 <- lapply(files2, function(f) readRDS(f)) 

test<-results_list_case5[[1]]

get_mode <- function(x) {
  uniq_x <- unique(x)
  uniq_x[which.max(tabulate(match(x, uniq_x)))]
}

summarize_vicatmix_output<-function(fits, true_alloc){
  map_dfr(seq_along(fits), function(r) {
    fit <- fits[[r]]
    z_true <- true_alloc[[r]]
    z_hat  <- fit$classification 
    K_hat<-get_mode(fit$clustersnum)
    
    tibble(rep = r,
           ARI = adjustedRandIndex(z_true, z_hat),
           K_hat = K_hat,
           PIP= I(list(fit$pip)))
  })
}


res_case5 <- summarize_vicatmix_output(
  fits = results_list_case5,
  true_alloc = alloclist5)


res_case7 <- summarize_vicatmix_output(
  fits = results_list_case7,
  true_alloc = alloclist7)

res_case5$scen<- "case5"
res_case7$scen<-"case7"

res_all <- bind_rows(res_case5, res_case7)
res_all$K_true<-5
res_all$method<- "VICatMix"

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

#Seems like we may need more shrinkage in the mixing weights prior



#Now we explore PIP output
extract_pip <- function(fit, rep, scen){
  tibble(
    rep = rep,
    scen = scen,
    var = paste0("X", seq_along(fit$pip)),
    pip = fit$pip
  )
}


pip_df <- bind_rows(
  map2_dfr(results_list_case5, seq_along(results_list_case5),
           ~ extract_pip(.x, .y, "case5")),
  map2_dfr(results_list_case7, seq_along(results_list_case7),
           ~ extract_pip(.x, .y, "case7"))
)



true_vars <- paste0("X", 1:12)
noise_vars <- paste0("X", 13:20)

pip_summary <- pip_df %>%
  group_by(scen, var) %>%
  summarise(
    mean_pip = mean(pip),
    sd_pip = sd(pip),
    .groups = "drop"
  ) %>%
  mutate(
    type = case_when(
      var %in% true_vars ~ "signal",
      TRUE ~ "noise"
    )
  )
ggplot(filter(pip_summary, scen == "case7"),
       aes(x = reorder(var, mean_pip), y = mean_pip, fill = type)) +
  geom_col() +
  coord_flip() +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme_minimal()



#Save results to home folder
saveRDS(res_all,file = "/n/home03/crodriguezcabrera/Mbeta_Project/Sim_results/vitcatmix_sim_tertiles.rds")




#VitCatmixAvg
files_path<-"/n/netscratch/stephenson_lab/Lab/crodriguez/OMbeta_outputs/"
files1 <- list.files(files_path, full.names = TRUE, pattern = "^VICatMixAvg_case5_rep.*\\.rds$")
files2 <- list.files(files_path, full.names = TRUE, pattern = "^VICatMixAvg_case7_rep.*\\.rds$")


results_list_case5 <- lapply(files1, function(f) readRDS(f)) 
results_list_case7 <- lapply(files2, function(f) readRDS(f)) 

test<-results_list_case5[[1]]
table(test$classification)


summarize_vicatmixAvg_output<-function(fits, true_alloc){
  map_dfr(seq_along(fits), function(r) {
    fit <- fits[[r]]
    z_true <- true_alloc[[r]]
    z_hat  <- fit$classification 
    K_hat<- max(fit$classification)
    
    tibble(rep = r,
           ARI = adjustedRandIndex(z_true, z_hat),
           K_hat = K_hat,
           PIP= I(list(fit$pip)))
  })
}


res_case5 <- summarize_vicatmixAvg_output(
  fits = results_list_case5,
  true_alloc = alloclist5)

res_case7 <- summarize_vicatmixAvg_output(
  fits = results_list_case7,
  true_alloc = alloclist7)

res_case5$scen<- "case5"
res_case7$scen<-"case7"

res_all <- bind_rows(res_case5, res_case7)
res_all$K_true<-5
res_all$method<- "VICatMixAvg"

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



#Save results to home folder
saveRDS(res_all,file = "/n/home03/crodriguezcabrera/Mbeta_Project/Sim_results/vitcatmixAvg_sim.rds")












