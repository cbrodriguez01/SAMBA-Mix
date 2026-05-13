library(tidyverse)
library(dplyr)
#Simulation results: comparative models combined
setwd("~/Mbeta_Project/Sim_results")

#Load all sim results
mclust_res<- readRDS("mclust_sim.rds")
clustvarsel_res<-readRDS("clustvarsel_sim.rds")
bayesbinmix_res<-readRDS("bayesbinmix_sim.rds")
vitcatmix_res<-readRDS("vitcatmix_sim_tertiles.rds")
sambo_acs<-readRDS("SAMBO_acsmock_case5_7_3.20.26.rds")#res_all and pip_df
sambo_res<-sambo_acs[[1]]
sambo_res$method<- "SAMBA-Mix"

#rename to match across
str(mclust_res)
mclust_res<- mclust_res %>% rename(scen = case)
clustvarsel_res<- clustvarsel_res %>% rename(scen = case)
sambo_res<-sambo_res %>%  rename( K_hat = Kmode)



#Keep: scen, rep, ARI,K_hat, method
keep<- c("scen", "rep", "ARI", "K_hat", "method")
mclust_res1<- mclust_res %>%  select(all_of(keep))
clustvarsel_res1<- clustvarsel_res %>%  select(all_of(keep))
bayesbinmix_res1<- bayesbinmix_res %>%  select(all_of(keep))
vitcatmix_res1<- vitcatmix_res %>%  select(all_of(keep))
#vitcatmixAvg_res1<- vitcatmixAvg_res %>%  select(all_of(keep))
sambo_res1<- sambo_res %>% select(all_of(keep))
#Combine
res_all<- bind_rows(mclust_res1, clustvarsel_res1, bayesbinmix_res1, vitcatmix_res1, sambo_res1)





#Plots for ARI and Khat
res_all$method <- factor(
  res_all$method,
  levels = c("BayesBinMix", "VICatMix", "mclust", "clustvarsel", "SAMBA-Mix")
)

method_cols <- c(
  "BayesBinMix" = "#D55E00",  # muted vermillion
  "clustvarsel" = "#009E73",  # bluish green
  "mclust"      = "#0072B2",  # blue
  "VICatMix"   = "#CC7",  # purple
  "SAMBA-Mix" = "#CC79A7"
)

# ggplot(res_all, aes(x=method, y=K_hat, fill=method)) +
#   geom_boxplot() +
#   facet_grid(~scen) +
#   theme_minimal()


P1<-ggplot(res_all, aes(x = method, y = K_hat, fill = method)) +
  geom_boxplot(
    width = 0.6,
    outlier.size = 0.8,
    outlier.alpha = 0.5,
    linewidth = 0.4
  ) +
  facet_wrap(~ scen, nrow = 1,
             labeller = labeller(scen = c(
               case5 = "Scenario 5",
               case7 = "Scenario 6"
             ))) +
  geom_hline(yintercept = 5,
             linetype = "dashed",
             linewidth = 0.4,
             color = "grey40") +
  scale_fill_manual(values = method_cols) +
  labs(
    x = NULL,
    y = expression(hat(K))) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    axis.line = element_line(linewidth = 0.4),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )



P2<-ggplot(res_all, aes(x = method, y = ARI, fill = method)) +
  geom_boxplot(
    width = 0.6,
    outlier.size = 0.8,
    outlier.alpha = 0.5,
    linewidth = 0.4
  ) +
  facet_wrap(~ scen, nrow = 1,
             labeller = labeller(scen = c(
               case5 = "Scenario 5",
               case7 = "Scenario 6"
             )))  +
  scale_fill_manual(values = method_cols) +
  labs(
    x = NULL,
    y = "Adjusted Rand Index (ARI)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    axis.line = element_line(linewidth = 0.4),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


#make a grid instead
#res_all_long<- pivot_longer(res_all, cols = c("ARI", "K_hat"), names_to = "measure", values_to = "Estimate")
# 
# ggplot(res_all,  aes(x = method, y = ARI, fill = method)) +
#   geom_boxplot(
#     width = 0.6,
#     outlier.size = 0.8,
#     outlier.alpha = 0.5,
#     linewidth = 0.4
#   ) +
#   facet_grid(scen ~ method, 
#              labeller = labeller(scen = c(
#                case5 = "Mock ACS case 1",
#                case7 = "Mock ACS case 2"
#              )))  +
#   scale_fill_manual(values = method_cols) +
#   labs(
#     x = NULL,
#     y = "Adjusted Rand Index (ARI)"
#   ) +
#   theme_classic(base_size = 12) +
#   theme(
#     legend.position = "top",
#     legend.title = element_blank(),
#     strip.background = element_rect(fill = "grey95", color = NA),
#     strip.text = element_text(face = "bold"),
#     axis.line = element_line(linewidth = 0.4)
#   )
# 


png(filename = "~/bayesmbmm/simulations/figures/Khat_allmethods2_SAMBAm.png")
P1
dev.off()

png(filename = "~/bayesmbmm/simulations/figures/ARI_allmethods2SAMBAm.png")
P2
dev.off()


#Models that do variable selection: VitCatMix and clustvarsel

##Filter


reps_withK5<-c(11, 22, 30, 34, 36, 39,4, 40, 45, 46,
               47, 52, 60, 62, 64, 70, 73,74,
               76, 79, 80, 87,9, 93, 95, 17, 22, 27 , 3, 41, 74, 78, 85)

res_all_5<- res_all %>% filter(rep %in% reps_withK5 & scen == "case5")
res_all_7<- res_all %>% filter(rep %in% reps_withK5 & scen == "case7")

ggplot(res_all_5, aes(x = method, y = K_hat, fill = method)) +
  geom_boxplot(
    width = 0.6,
    outlier.size = 0.8,
    outlier.alpha = 0.5,
    linewidth = 0.4
  ) +
  facet_wrap(~ scen, nrow = 1,
             labeller = labeller(scen = c(
               case5 = "Scenario 5",
               case7 = "Scenario 6"
             ))) +
  geom_hline(yintercept = 5,
             linetype = "dashed",
             linewidth = 0.4,
             color = "grey40") +
  scale_fill_manual(values = method_cols) +
  labs(
    x = NULL,
    y = expression(hat(K))) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    axis.line = element_line(linewidth = 0.4),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )









