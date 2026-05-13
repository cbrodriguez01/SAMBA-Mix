#"AHRQ dataset (census tract 2020)"
#Last I updated this code: 2/7/2026


#PENDING
#  Is  this data is using 2010 or 2020 Census boundaries, as the numbers changed significantly?
# Consider imputation strategies for missing data in the future.

#Load packages
library(readxl)
library(tidyverse)
library(ggplot2)
library(GGally)
library(visdat)
library(psych)
library(corrplot)
library(gridExtra)
#library(clustvarsel)
library(patchwork)

# Key variables
# 
# -   YEAR: SDOH file year
# -   TRACTFIPS: State-county-census tract FIPS Code (11-digit)
# -   COUNTYFIPS: State-county FIPS Code (5-digit)
# -   STATEFIPS: State FIPS Code (2-digit)
# -   STATE: State name
# -   COUNTY: County name
# -   REGION: Census region name
# -   TERRITORY: Territory indicator (1= U.S. Territory, 0= U.S. State or DC)
# 


sdoh_ahrq<-read_xlsx("~/AHRQData/sdoh_2020_tract_1_0.xlsx", sheet = "Data")
varlist <- c(
  "TRACTFIPS", "COUNTYFIPS", "STATEFIPS", "STATE", "COUNTY", "REGION",
  "ACS_PCT_FOREIGN_BORN",
  "ACS_PCT_ENGL_NOT_ALL",
  "ACS_PCT_HH_LIMIT_ENGLISH",
  "ACS_PCT_HH_NO_COMP_DEV",
  "ACS_PCT_HH_NO_INTERNET",
  "ACS_PCT_CHILD_1FAM",
  "ACS_PCT_UNEMPLOY",
  "ACS_PCT_INC50",
  "ACS_PCT_HH_PUB_ASSIST",
  "ACS_PCT_LT_HS",
  "ACS_PCT_1UP_RENT_1ROOM",
  "ACS_PCT_RENTER_HU",
  "ACS_PCT_HU_NO_VEH",
  "ACS_PCT_MEDICAID_ANY",
  "ACS_PCT_GRP_QRT",
  "ACS_PCT_VET",
  "ACS_PCT_HH_SMARTPHONE_ONLY", 
  "ACS_PCT_DISABLE"
)

#FIPS: 25 & 36
sdoh_ahrq_statesCTs<- sdoh_ahrq %>% filter(TERRITORY == 0) %>% 
  filter(STATE %in% c("Massachusetts", "New York")) %>% select(all_of(varlist))


vis_value(sdoh_ahrq_statesCTs %>% select(-c("TRACTFIPS", "COUNTYFIPS", "STATEFIPS", "STATE", "COUNTY", "REGION")))

#Missingness for MA
vis_miss(sdoh_ahrq_statesCTs %>% filter(STATEFIPS == "25"), warn_large_data = FALSE, sort_miss = TRUE )
#Missingness for NY
vis_miss(sdoh_ahrq_statesCTs %>% filter(STATEFIPS == "36"), warn_large_data = FALSE, sort_miss = TRUE)

#NICE! it's very low, so we can just stick to complete cases for this analysis


sdoh_ahrq_statesCTs_complete <- sdoh_ahrq_statesCTs %>% drop_na()


##------------ Overall distribution of the variables --------------------------

numeric_vars <- sdoh_ahrq_statesCTs_complete %>% select(where(is.numeric)) %>% names()

sdoh_ahrq_states_long<- sdoh_ahrq_statesCTs_complete %>% pivot_longer(7:28, names_to = "SDoH_variable",values_to = "value")
selected_vars <- c(
  "ACS_PCT_FOREIGN_BORN",
  "ACS_PCT_ENGL_NOT_ALL",
  "ACS_PCT_HH_NO_INTERNET",
  "ACS_PCT_UNEMPLOY",
  "ACS_PCT_INC50",
  "ACS_PCT_LT_HS",
  "ACS_PCT_RENTER_HU",
  "ACS_PCT_MEDICAID_ANY"
)

sdoh_ahrq_states_long_selected <- sdoh_ahrq_statesCTs_complete %>%
  pivot_longer(cols = all_of(selected_vars), names_to = "SDoH_variable", values_to = "value")

p1<-sdoh_ahrq_states_long_selected %>% filter(STATEFIPS == 36) %>% ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.7) +
  facet_wrap(~SDoH_variable, scales = "free") +
  theme_minimal() + theme(
    plot.title = element_text(size = 17),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 8)) +
  labs(title = "Distribution of Selected SDoH Variables 2020 (NY)",
       x = "Proportions", y = "Count")


p2<-sdoh_ahrq_states_long_selected %>% drop_na() %>% filter(STATEFIPS == 25) %>% ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.7) +
  facet_wrap(~SDoH_variable, scales = "free") +
  theme_minimal() + theme(
    plot.title = element_text(size = 17),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 8)) +
  labs(title = "Distribution of Selected SDoH Variables 2020 (MA)",
       x = "Proportions", y = "Count")

v<- c(
  "ACS_PCT_FOREIGN_BORN",
  "ACS_PCT_ENGL_NOT_ALL",
  "ACS_PCT_UNEMPLOY",
  "ACS_PCT_LT_HS",
  "ACS_PCT_RENTER_HU",
  "ACS_PCT_MEDICAID_ANY"
)
#https://stackoverflow.com/questions/68553330/r-ggpairs-change-the-colour-of-bar-histograms-on-diagonal
#MA
MAdat<-sdoh_ahrq_statesCTs %>% filter(STATEFIPS == 25) 
cor_matMA <- cor(MAdat[numeric_vars], use = "pairwise.complete.obs")
corrplot(cor_matMA, method = "color", type = "upper", tl.cex = 0.6, title = "Correlation Matrix for MA SDoH Variables", mar = c(0,0,1,0))


pp_MA<-MAdat %>% select(v) %>% ggpairs(
  upper = list(continuous = wrap("cor", size = 3.5)),
  lower = list(continuous = wrap("points", alpha = 0.3, size = 0.5)),
  diag = list(continuous = wrap("barDiag", bins = 30, fill = "steelblue", color = "white", alpha = 0.7))) + 
  theme(strip.text = element_text(size = 4)) 


#NY
NYdat<-sdoh_ahrq_statesCTs %>% filter(STATEFIPS == 36)
cor_matNY <- cor(NYdat[numeric_vars], use = "pairwise.complete.obs")
corrplot(cor_matNY, method = "color", type = "upper", tl.cex = 0.6, title = "Correlation Matrix for NY SDoH Variables", mar = c(0,0,1,0))

pp_NY<-NYdat %>% select(v) %>% ggpairs(
  upper = list(continuous = wrap("cor", size = 3.5)),
  lower = list(continuous = wrap("points", alpha = 0.3, size = 0.5)),
  diag = list(continuous = wrap("barDiag", bins = 30, fill = "steelblue", color = "white", alpha = 0.7))) + 
  theme(strip.text = element_text(size = 4)) 


png("figures/pairwise_plots_MA.png", width = 1000, height = 1000, res = 150)
pp_MA
dev.off()

png("figures/pairwise_plots_NY.png", width = 1000, height = 1000, res = 150)
pp_NY
dev.off()


#Boxplots to compare distributions between MA and NY for a subset of variables

boxplot_comparison <- sdoh_ahrq_states_long_selected %>% drop_na() %>%
  ggplot(aes(x = STATE, y = value, fill = STATE)) +
  geom_boxplot() +
  facet_wrap(~SDoH_variable, scales = "free") +
  theme_minimal() + theme(
    plot.title = element_text(size = 17),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 10),
    legend.position = "none"
  ) +
  labs(title = "Comparison of Selected SDoH Variables between MA and NY (2020)",
       x = "State", y = "Proportions")



pdf("SDoH_Distributions_MA_NY.pdf", width = 14, height = 10)
p1
p2
#boxplot_comparison
dev.off()


#Median and  IQRs
summary_stats <- sdoh_ahrq_statesCTs %>%
  group_by(STATE) %>%
  summarise(across(all_of(numeric_vars), list(median = ~median(.x, na.rm = TRUE), IQR = ~IQR(.x, na.rm = TRUE)), .names = "{col}_{fn}"))

#Output as CSV file
write.csv(summary_stats, "summary_stats_MA_NY.csv", row.names = FALSE)
#use table1 package later for better formatting



#Save data as Rdata
#convert percentages into decimals for easier interpretation in the future
sdoh_ahrq_statesCTs_complete1 <- sdoh_ahrq_statesCTs_complete %>%
  mutate(across(all_of(numeric_vars), ~ .x / 100))

NYdat<-sdoh_ahrq_statesCTs_complete1 %>% filter(STATEFIPS == 36) %>% select(-c(1:6))
MAdat<-sdoh_ahrq_statesCTs_complete1 %>% filter(STATEFIPS == 25) %>% select(-c(1:6))

eps <- 1e-3
NYdat_adj <- as.data.frame(
  lapply(NYdat, function(x) {
    x * (1 - 2 * eps) + eps
  })
)

MAdat_adj <- as.data.frame(
  lapply(MAdat, function(x) {
    x * (1 - 2 * eps) + eps
  })
)


saveRDS(list(NYdat_adj, MAdat_adj), file = "~/AHRQData/sdoh_ahrq_statesCTs_complete.rds")


#save with cts
MAdatcts<-sdoh_ahrq_statesCTs_complete1 %>% filter(STATEFIPS == 25) %>% select(c(1:6))
NYdatcts<-sdoh_ahrq_statesCTs_complete1 %>% filter(STATEFIPS == 36)  %>% select(c(1:6))


saveRDS(list(NYdatcts, MAdatcts), file = "~/AHRQData/sdoh_ahrq_statesCTs.rds")



#A few notes
#  many variables where:
#   Min = 0.001
# 1st quartile also near 0.001
# Examples:
#   ACS_PCT_ENGL_NOT_ALL
# ACS_PCT_1UP_RENT_1ROOM
# ACS_PCT_GRP_QRT
#   These variables are quasi-zero-inflated




