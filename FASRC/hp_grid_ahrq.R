library(dplyr)

e0_vec <- c(0.005)
rho_mean <- c(0.85, 0.6)
conc <- 10

hp <- tibble(
  rho_mean = rho_mean,
  a_rho = rho_mean * conc,
  b_rho = (1 - rho_mean) * conc
)

hp_all <- expand.grid(
  e0 = e0_vec,
  rho_mean = rho_mean,
  stringsAsFactors = FALSE
) %>%
  left_join(hp, by = "rho_mean")

saveRDS(hp_all, "~/AHRQData/hp_grid.rds")

write.csv(hp_all, "~/AHRQData/hp_grid.csv")
