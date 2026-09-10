#Going over model output-- updated in September 9, 2026 to add map for all census tracts
setwd("~/bayesmbmm/results")
library(nimble)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(coda)
library(bayesplot)
library(flextable)
library(officer)
library(writexl)
library(tigris) # To get MA shapefiles
library(sf)
library(ggtext)

ahrq<-readRDS(file = "~/AHRQData/sdoh_ahrq_statesCTs_complete.rds")
madat<-ahrq[[2]]
features<-colnames(madat)
#make names nicer
feature_names <- c(
  "Foreign Born", "English Not All", "Limited English", 
  "No Comp Dev", "No Internet", "Single-parent", 
  "Unemployment", "Low Income", "Public Assistance", 
  "Low Ed", "Rent-overcrowd", "Renter-occupied", 
  "No Vehicle", "Medicaid", "Group Qrt", 
  "Veteran", "Smartphone Only", "Disabled"
)

#res <- readRDS("/n/netscratch/stephenson_lab/Lab/crodriguez/AHRQ_application/ahrqMA_sambo_3.27.26.rds")
res <- readRDS("/n/netscratch/stephenson_lab/Lab/crodriguez/AHRQ_application/ahrqMA_sambamix_8.28.26.rds")

fit <- res[[1]]
runtime <- res[[2]]/3600 # 10 hours
Kmode<-fit$K

active_feats <- which(fit$FeatureSal_out[[2]] > 0.5)
p<-ncol(madat)
# 
mcmc_samples <- fit$modeloutput_perchain$parameters.ECR.mcmc
#Based on the trace plots, we are going to add more burn in and remove ~ 2000 
#samples more so that we have niter=20k and burnin=9k

mcmc_samples_burned<-mcmc_samples[-c(1:3000), ]
saveRDS(mcmc_samples_burned, file = "~/Mbeta_Project/mcmc_samples_burned_SAMBAmix_AHRQ_8.28.26.rds")


sum_stats<-samplesSummary(mcmc_samples_burned)

# lambda[1]  17.62618913  17.62154351 0.3116831706  17.01547864  18.21665968
# lambda[2]  11.59738675  11.59532723 0.2323263330  11.14913318  12.04890148
# lambda[3]   2.91663935   2.91760514 0.1086787033   2.70512502   3.12052438


mu_rows <- grep("^mu\\[", rownames(sum_stats))
mu_summary <- sum_stats[mu_rows, ]

# 3. Reshape from a long list into a 3x18 matrix (3 clusters, 18 variables)
# mu[1,1], mu[2,1], mu[3,1], mu[1,2]... 
# The loop ensures we map the right 'mu' to the right cluster/feature

mu_matrix <- matrix(NA, nrow = Kmode, ncol = p)

for(k in 1:Kmode) {
  for(j in 1:p) {
    row_name <- paste0("mu[", k, ", ", j, "]")
    if(row_name %in% rownames(mu_summary)) {
      mu_matrix[k, j] <- mu_summary[row_name, "Mean"]
    }
  }
}

colnames(mu_matrix) <- feature_names


# Convert your matrix to a dataframe
cluster_tables <- list()
for(k in 1:3) {
  # Identify the rows for this specific cluster's mu parameters
  # Format in sum_stats is "mu[k, 1]", "mu[k, 2]", etc.
  rows <- paste0("mu[", k, ", ", 1:18, "]")
  # Extract and rename
  tab <- as.data.frame(sum_stats[rows, c("Mean", "St.Dev.", "95%CI_low", "95%CI_upp")])
  colnames(tab) <- paste0(c("Mean", "SD", "LCI", "UCI"), "_C", k)
  
  cluster_tables[[k]] <- tab
}

# Combine into one wide table and add variable names
full_table <- do.call(cbind, cluster_tables)
rownames(full_table) <- feature_names


#i want to get word doc table

# Function to extract and format: "Mean (LCI-UCI)"
get_formatted_cluster <- function(k, stats, names) {
  rows <- paste0("mu[", k, ", ", 1:18, "]")
  m <- stats[rows, "Mean"]
  l <- stats[rows, "95%CI_low"]
  u <- stats[rows, "95%CI_upp"]
  
  paste0(sprintf("%.3f", m), " (", sprintf("%.4f", l), "-", sprintf("%.4f", u), ")")
}

word_df <- data.frame(
  Variable = feature_names,
  Cluster_1 = get_formatted_cluster(1, sum_stats, feature_names),
  Cluster_2 = get_formatted_cluster(2, sum_stats, feature_names),
  Cluster_3 = get_formatted_cluster(3, sum_stats, feature_names),
  stringsAsFactors = FALSE
)

# Create the Flextable
ft <- flextable(word_df) %>%
  set_header_labels(
    Variable = "SDOH Variable",
    Cluster_1 = "Cluster 1(n=885)",
    Cluster_2 = "Cluster 2(n=516)",
    Cluster_3 = "Cluster 3(n=177)"
  ) %>%
  autofit() %>%
  # Add a title
  add_header_lines("Table 1: Posterior Mean SDOH Profiles with 95% Credible Intervals") %>%
  # Visual formatting
  bold(part = "header") %>%
  bg(j = 3, bg = "#F2F2F2", part = "body") %>% # Highlight the High-Need Cluster
  border_remove() %>%
  hline_top(part = "all", border = fp_border(width = 2)) %>%
  hline_bottom(part = "body", border = fp_border(width = 2))

#  Export to Word
read_docx() %>%
  body_add_flextable(ft) %>%
  print(target = "~/bayesmbmm/results/figures_tabs/AHRQ_MA_Cluster_Results_9.9.26.docx")



#We want to group by domain to be able to plot
var_info <- data.frame(
  Variable = c(
    # Economic
    "Unemployment", "Low Income", "Public Assistance",
    
    # Education
    "Low Ed",
    
    # Healthcare
    "Medicaid",
    
    # Physical infrastructure
    "No Comp Dev", "No Internet", "Rent-overcrowd",
    "Renter-occupied", "No Vehicle", "Group Qrt", "Smartphone Only",
    
    # Social
    "Foreign Born", "English Not All", "Limited English",
    "Single-parent", "Veteran", "Disabled"
  ),
  Domain = c(
    rep("Economic context", 3),
    "Education",
    "Healthcare context",
    rep("Physical infrastructure", 7),
    rep("Social context", 6)
  )
)



# Extract Mean, LCI, UCI for all clusters
plot_data <- data.frame()
for(k in 1:3) {
  rows <- paste0("mu[", k, ", ", 1:18, "]")
  temp <- as.data.frame(sum_stats[rows, c("Mean", "95%CI_low", "95%CI_upp")])
  temp$Variable <- feature_names
  temp$Cluster <- paste("Profile", k)
  temp$PIP <- as.numeric(fit$FeatureSal_out[[2]])
  plot_data <- rbind(plot_data, temp)
}

colnames(plot_data)[2:3] <- c("LCI", "UCI")
# Filter to 'Saliant' features (PIP > 0.5) for cleaner plots
plot_data_active <- plot_data %>% filter(PIP > 0.5)

plot_data_active <- plot_data_active %>%
  left_join(var_info, by = "Variable")

# enforce exact table order
plot_data_active$Variable <- factor(
  plot_data_active$Variable,
  levels = var_info$Variable
)

#forest plot

fp <- ggplot(plot_data_active, aes(x = Mean, y = Variable, color = Cluster,shape = Cluster)) +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI),
                 height = 0.3, size = 0.8,
                 position = position_dodge(0.6)) +
  geom_point(size = 2.8, position =  position_dodge(width = 0.6)
) +
  facet_grid(Domain ~ ., scales = "free_y", space = "free_y") +

  scale_color_manual(values = c("#2c7bb6", "#d7191c", "#fdae61")) +
  scale_shape_manual(values = c(16, 17, 18)) +

  theme_minimal() +
  theme(
    strip.text.y = element_text(angle = 0, face = "bold"),
    panel.spacing = unit(0.8, "lines")
  ) +

  labs(x = "Proportion", y = "")
  
  
  
  
  #theme_minimal() +
  #scale_color_manual(values = c("#2c7bb6", "#d7191c", "#fdae61")) +
  #labs(
    #title = "Posterior Means  95% Credible Intervals",
       #subtitle = "Non-overlapping bars indicate statistically distinct SDOH profiles",
      # x = "Proportion", y = "")

ggsave("figures_tabs/ForestPlot_MA_AHRQ.pdf", fp, width = 8, height = 6)
#title = "Posterior Means  95% Credible Intervals",
#subtitle = "Non-overlapping bars indicate statistically distinct SDOH profiles",

# #radar plot
# rp <- ggplot(plot_data_active, aes(x = Variable, y = Mean, group = Cluster, color = Cluster)) +
#   geom_polygon(fill = NA, size = 1) +
#   coord_polar() +
#   facet_wrap(~Cluster) +
#   theme_bw() +
#   theme(axis.text.x = element_text(size = 7)) +
#   labs(title = "Community SDOH Signatures", y = "")
# 
# ggsave("RadarPlot_AHRQ.pdf", rp, width = 10, height = 5)
# 

#I want to add more structure-- and move domains to the left side
y_levels <- c(
  "Economic context",
  "Unemployment", "Low Income", "Public Assistance",
  
  "Education",
  "Low Ed",
  
  "Healthcare context",
  "Medicaid",
  
  "Physical infrastructure",
  "No Comp Dev", "No Internet", "Rent-overcrowd",
  "Renter-occupied", "No Vehicle", "Group Qrt", "Smartphone Only",
  
  "Social context",
  "Foreign Born", "English Not All", "Limited English",
  "Single-parent", "Veteran", "Disabled"
)

headers <- c(
  "Economic context", "Education", "Healthcare context",
  "Physical infrastructure", "Social context"
)

headers_df <- data.frame(
  Variable = headers,
  Mean = NA, LCI = NA, UCI = NA,
  Cluster = NA
)

plot_data_final <- bind_rows(plot_data_active, headers_df)

plot_data_final$Variable <- factor(
  plot_data_final$Variable,
  levels = rev(y_levels)
)



plot_data_final <- plot_data_final %>%
  mutate(
    is_header = Variable %in% headers,
    label = Variable
  )

p_left <- ggplot(plot_data_final, aes(y = Variable)) +
  geom_text(aes( x = ifelse(is_header,0,0.25),
      label = Variable,
      fontface = ifelse(is_header, "bold", "plain")
    ),
    hjust = 0,
    size = 3.8,
    lineheight = 0.2
  ) +
  scale_y_discrete(
    limits = levels(plot_data_final$Variable),
    drop = FALSE
  ) +
  xlim(0, 1) +   # ← IMPORTANT: gives text room to render
  theme_void() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
  
  
  #theme(plot.margin = margin(5, 5, 5, 5))


p_right <- ggplot(plot_data_final, aes(x = Mean, y = Variable, color = Cluster, shape = Cluster)) +
  geom_errorbarh(
    data = subset(plot_data_final, !is.na(Mean)),
    aes(xmin = LCI, xmax = UCI),
    height = 0.3,
    size = 0.5,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    data = subset(plot_data_final, !is.na(Mean)),
    size = 2.8,
    position = position_dodge(width = 0.6)
  ) +
  scale_y_discrete(limits = levels(plot_data_final$Variable)) +
  scale_color_manual(name = "NSDoH Profile", values = c("#2c7bb6", "#d7191c", "#fdae61")) +
  scale_shape_manual(name = "NSDoH Profile", values = c(16, 17, 18), na.translate = FALSE) +
  theme_minimal(base_size = 18) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(5, 5, 5, 0),
  ) + labs(x = "Mean proportion")

final_plot <- p_left + p_right + plot_layout(widths = c(2, 5))


ggsave("figures_tabs/ForestPlot_MA_AHRQ_domain.pdf", final_plot)



########different format##########

# ---------------------------------------------------------
# Y-axis structure
# ---------------------------------------------------------

y_levels <- c(
  "Economic context",
  "Unemployment", "Low Income", "Public Assistance",
  
  "Education",
  "Low Ed",
  
  "Healthcare context",
  "Medicaid",
  
  "Physical infrastructure",
  "No Comp Dev", "No Internet", "Rent-overcrowd",
  "Renter-occupied", "No Vehicle", "Group Qrt", "Smartphone Only",
  
  "Social context",
  "Foreign Born", "English Not All", "Limited English",
  "Single-parent", "Veteran", "Disabled"
)

headers <- c(
  "Economic context",
  "Education",
  "Healthcare context",
  "Physical infrastructure",
  "Social context"
)

# ---------------------------------------------------------
# Add blank rows for domain headings
# ---------------------------------------------------------

headers_df <- data.frame(
  Variable = headers,
  Mean = NA,
  LCI = NA,
  UCI = NA,
  Cluster = NA
)

plot_data_final <- bind_rows(
  plot_data_active,
  headers_df
)

plot_data_final$Variable <- factor(
  plot_data_final$Variable,
  levels = rev(y_levels)
)

plot_data_final <- plot_data_final %>%
  mutate(
    is_header = Variable %in% headers,
    label = Variable
  )


# ---------------------------------------------------------
# LEFT PANEL
# Domain headings + variable labels
# ---------------------------------------------------------

p_left <- ggplot(
  plot_data_final,
  aes(y = Variable)
) +
  geom_text(
    aes(
      x = ifelse(is_header, 0, 0.25),
      label = Variable,
      fontface = ifelse(is_header, "bold", "plain")
    ),
    hjust = 0,
    size = 3.8
  ) +
  scale_y_discrete(
    limits = levels(plot_data_final$Variable),
    drop = FALSE
  ) +
  xlim(0, 1) +
  theme_void() +
  theme(
    plot.margin = margin(
      t = 5,
      r = 0,
      b = 5,
      l = 5
    )
  )


# ---------------------------------------------------------
# RIGHT PANEL
# One panel per profile
# ---------------------------------------------------------

p_right <- ggplot(
  data = plot_data_final %>% filter(!is.na(Mean)),
  aes(
    x = Mean,
    y = Variable,
    color = Cluster,
    shape = Cluster
  )
) +
  
  # Confidence intervals
  geom_errorbarh(
    aes(
      xmin = LCI,
      xmax = UCI
    ),
    height = 0.25,
    linewidth = 0.7
  ) +
  
  # Posterior means
  geom_point(
    size = 2
  ) +
  
  # Separate panel for each profile
  facet_grid(
    . ~ Cluster,
    scales = "free_x"
  ) +
  
  scale_y_discrete(
    limits = levels(plot_data_final$Variable),
    drop = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      "#2c7bb6",
      "#d7191c",
      "#fdae61"
    ),
    guide = "none"
  ) +
  
  scale_shape_manual(
    values = c(
      16,
      17,
      18
    ),
    guide = "none"
  ) +
  
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  
  labs(
    x = "Mean proportion",
    y = NULL
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    # remove duplicated y labels
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    # facet titles
    strip.text = element_text(
      size = 14,
      face = "bold"
    ),
    
    strip.background = element_rect(
      fill = "grey92",
      color = "grey70",
      linewidth = 0.5
    ),
    
    # emphasize horizontal structure
    panel.grid.major.y = element_line(
      color = "grey90",
      linewidth = 0.5
    ),
    
    panel.grid.minor = element_blank(),
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 0
    )
  )


# ---------------------------------------------------------
# COMBINE LABEL PANEL + PROFILE PANELS
# ---------------------------------------------------------

final_plot <- p_left + p_right +
  plot_layout(
    widths = c(2.3, 7)
  )

final_plot





ggsave("figures_tabs/ForestPlot_MA_AHRQ_domain_cols.pdf", final_plot)


##########



#lollipop
# Calculate 'Impact' as the difference between Cluster 2 and the State average (Cluster 1)
c1c2_data <- plot_data_active %>%
  filter(Cluster %in% c("Cluster 1", "Cluster 2")) %>%
  select(Variable, Cluster, Mean) %>%
  pivot_wider(names_from = Cluster, values_from = Mean) %>%
  mutate(Diff = `Cluster 2` - `Cluster 1`)

c1c3_data <- plot_data_active %>%
  filter(Cluster %in% c("Cluster 1", "Cluster 3")) %>%
  select(Variable, Cluster, Mean) %>%
  pivot_wider(names_from = Cluster, values_from = Mean) %>%
  mutate(Diff = `Cluster 3` - `Cluster 1`)

lp1 <- ggplot(c1c2_data, aes(x = reorder(Variable, Diff), y = Diff)) +
  geom_segment(aes(xend = Variable, yend = 0), color = "grey") +
  geom_point(size = 4, color = "#d7191c") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top drivers of  Cluster 2",
       subtitle = "Difference in Proportion(vs. Cluster 1)",
       x = "", y = "Increase in Proportion")

lp2 <- ggplot(c1c3_data, aes(x = reorder(Variable, Diff), y = Diff)) +
  geom_segment(aes(xend = Variable, yend = 0), color = "grey") +
  geom_point(size = 4, color = "#d7191c") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top drivers of  Cluster 3",
       subtitle = "Difference in Proportion(vs. Cluster 1)",
       x = "", y = "Increase in Proportion")

ggsave("figures_tabs/LollipopPlotC2vsC1_MAAHRQ.pdf", lp1, width = 7, height = 5)
ggsave("figures_tabs/LollipopPlotC3vsC1_MAAHRQ.pdf", lp2, width = 7, height = 5)



#Export as csv file
write_xlsx(
  list(
    sheet1 = plot_data,
    sheet2 = c1c2_data,
    sheet3 = c1c3_data
  ),
  path = "figures_tabs/AHRQ_MA_results.xlsx"
)

#cluster assignment probabilities 
ecrclusters<-as.vector(fit$modeloutput_perchain$clusterMembershipPerMethod)
probsclust<- fit$modeloutput_perchain$classificationProbabilities.ECR
maxprob<-apply(probsclust, 1,max)
zdat<- as.data.frame(cbind(ecrclusters,maxprob))
zdat$clusters<- factor(zdat$ecrclusters,
                          levels = 1:3, labels = paste0("Cluster",sep= " ", 1:3))


box<-zdat %>%  ggplot(aes(x = clusters, y = maxprob, fill = clusters)) +
  geom_boxplot() +
  labs(fill = "Profile Assignment", 
       y = "Profile Assignment Probability", 
       x= "") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 16),
        #axis.title.x = element_text(size = 18, color = "black", face = "bold"),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size=12), 
        axis.text.x = element_text(size=12))



ggsave("figures_tabs/profassignprob_AHRQ_MA.pdf", box, width = 7, height = 5)



#################################################################################
#################################################################################
#################################################################################
#MAPPING
#MAP these profiles similar to paper 1
cts<-readRDS("~/AHRQData/sdoh_ahrq_statesCTs.rds")
madatcts<-cts[[2]]
macomp<-cbind(madatcts, madat)
macomp$cluster_assignment_ecr<-ecrclusters
#table(macomp$cluster_assignment_ecr, macomp$cluster_assignment_prob)#perfect

macomp<- macomp %>% rename(GEOID=TRACTFIPS)

#Download MA Shapefile (using tigris)
ma_shape <- tracts(state = "MA", cb = TRUE)
# 1. Filter Shapefile to the Boston Metro Counties
# 025=Suffolk (Boston), 017=Middlesex (Cambridge/Somerville), 
# 021=Norfolk (Quincy), 009=Essex (Lynn)

boston_metro <- ma_shape  %>% filter(COUNTYFP %in% c("025", "017", "021", "009"))
map_df <- boston_metro %>%
  right_join(macomp, by = "GEOID")
bmap<-ggplot(map_df) +
  geom_sf(aes(fill = as.factor(cluster_assignment_ecr)), color = NA) +
  scale_fill_manual(values = c("1" = "#2c7bb6", "2" = "#d7191c", "3" = "#fdae61")) + 
  theme_void() + 
  #labs(title = "Geographic Distribution of SDOH Clusters in Boston Metro Counties",
  #     subtitle = "Massachusetts Census Tracts (AHRQ 2020)")+ 
  annotate("text", x = -71.0589, y = 42.3601, label = "Boston", fontface = "bold") +
  annotate("text", x = -71.1097, y = 42.3736, label = "Cambridge", fontface = "italic") 
  #annotate("text", x = -71.07, y = 42.30, label = "Dorchester", size = 3)
ggsave("figures_tabs/bostonmetromap_MAAHRQ.pdf", bmap, width = 7, height = 5)

#UPDATED MAP: 9/9/26
map_df_all <- ma_shape %>%
  right_join(macomp, by = "GEOID")

statemap<-ggplot(map_df_all) +
  geom_sf(aes(fill = as.factor(cluster_assignment_ecr)), color = NA) +
  scale_fill_manual(values = c("1" = "#2c7bb6", "2" = "#d7191c", "3" = "#fdae61"),
                    labels = c("Profile 1", "Profile 2", "Profile 3"),
                    name = "NSDoH Profile") + theme_void() 
 # labs(title = "Geographic Distribution of NSDoH Profiles ",
  #     subtitle = "Massachusetts Census Tracts (AHRQ 2020)")
ggsave("figures_tabs/MAmapSept_MAAHRQ.pdf", statemap, width = 7, height = 5)

statemap + bmap +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


#################################################################################
#updated map-- with inset
#https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/

map_df_all <- ma_shape %>%
  right_join(macomp, by = "GEOID")

boston_xlim <- c(-71.35, -70.85)
boston_ylim <- c(42.15, 42.60)

locator_box <- data.frame(
  xmin = boston_xlim[1], xmax = boston_xlim[2],
  ymin = boston_ylim[1], ymax = boston_ylim[2]
)




statemap <- ggplot(map_df_all) +
  geom_sf(
    aes(fill = as.factor(cluster_assignment_ecr)),
    color = NA,
    linewidth = 0.05
  ) +
  # Locator box: shows the viewer exactly what area the inset zooms into
  geom_rect(
    data = locator_box,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA,
    color = "black",
    linewidth = 0.5,
    linetype = "dashed"
  ) +
  
  scale_fill_manual(
    values = c(
      "1" = "#2c7bb6",
      "2" = "#d7191c",
      "3" = "#fdae61"
    ),
    labels = c(
      "Profile 1",
      "Profile 2",
      "Profile 3"
    ),
    name = "NSDoH Profile"
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )


#Greater boston
bmap <- ggplot(map_df_all) +
  geom_sf(
    aes(fill = as.factor(cluster_assignment_ecr)),
    color = NA,
    linewidth = 0.08
  ) +
  scale_fill_manual(
    values = c(
      "1" = "#2c7bb6",
      "2" = "#d7191c",
      "3" = "#fdae61"
    )
  ) +
  
  # Zoom to Greater Boston
  coord_sf(
    xlim = c(-71.35, -70.85),
    ylim = c(42.15, 42.60),
    expand = FALSE
  ) +
  
  # annotate(
  #   "text",
  #   x = -71.0589,
  #   y = 42.3601,
  #   label = "Boston",
  #   fontface = "bold",
  #   size = 3
  # ) +
  # 
  # annotate(
  #   "text",
  #   x = -71.1097,
  #   y = 42.3736,
  #   label = "Cambridge",
  #   fontface = "italic",
  #   size = 3
  # ) +
  
  labs(
    title = "Greater Boston"
  ) +
  
  theme_void() +
  
  # Remove duplicate legend and add inset border
  theme(
    legend.position = "none",
    
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 3)
    ),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.7
    ),
    
    plot.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.7
    )
  )


# Add inset

final_map <- statemap +
  inset_element(
    bmap,
    
    # position of inset within statewide figure
    left   = 0.80,
    bottom = 0.49,
    right  = 0.98,
    top    = 0.97,
    align_to = "full"
  )


final_map

ggsave("figures_tabs/statewithinset_MAAHRQ.pdf", final_map, width = 7, height = 5)

#################################################################################
#################################################################################
#################################################################################
### Explore traceplots
# s <- as.matrix(fit$samples_raw)
# nk<-fit$nk_iter
# colMeans(nk) 
# active_clusts <- order(colMeans(fit$nk_iter), decreasing = TRUE)[1:Kmode]

## plot trace of pi
pi_samples<-mcmc_samples_burned[,paste0("pi[", 1:3, "]")]
color_scheme_set("brightblue")
mcmc_trace(pi_samples, 
           pars = paste0("pi[", 1:3, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for pi")

lambda_samples<-mcmc_samples[,paste0("lambda[", 1:3, "]")]
color_scheme_set("brightblue")
mcmc_trace(lambda_samples, 
           pars = paste0("lambda[", 1:3, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for lambda")


#MUCH BETTER after removing 3000 samps so that burnin=10k
lambda_samples1<-mcmc_samples_burned[,paste0("lambda[", 1:3, "]")]
mcmc_trace(lambda_samples1, 
           pars = paste0("lambda[", 1:3, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for lambda")


## plot trace of mu for feature

mu_samples1<-mcmc_samples_burned[,paste0("mu[", 1, ", ", 1:18, "]")]
color_scheme_set("brightblue")

mcmc_trace(mu_samples1, 
           pars = paste0("mu[", 1, ", ", 1:10, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for mu (cluster 1)")

mcmc_trace(mu_samples1, 
           pars = paste0("mu[", 1, ", ", 11:18, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for mu (cluster 1)")



mu_samples2<-mcmc_samples_burned[,paste0("mu[", 2, ", ", 1:18, "]")]
color_scheme_set("brightblue")
mcmc_trace(mu_samples2, 
           pars = paste0("mu[", 2, ", ", 1:10, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for mu (cluster 2)")

mcmc_trace(mu_samples2, 
           pars = paste0("mu[", 2, ", ", 11:18, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for mu (cluster 2)")



mu_samples3<-mcmc_samples_burned[,paste0("mu[", 3, ", ", 1:18, "]")]
color_scheme_set("brightblue")
mcmc_trace(mu_samples3, 
           pars = paste0("mu[", 3, ", ", 1:10, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for mu  (cluster 3)")
mcmc_trace(mu_samples3, 
           pars = paste0("mu[", 3, ", ", 11:18, "]"),
           facet_args = list(ncol = 2)) +
  ggtitle("Trace plots for mu  (cluster 3)")
