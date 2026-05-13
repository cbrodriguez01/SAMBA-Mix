library(nimble)
library(coda)

res <- readRDS("/n/netscratch/stephenson_lab/Lab/crodriguez/AHRQ_application/ahrqNY_sambo_4.8.26.rds")
fit <- res[[1]]
runtime <- res[[2]]/3600 # 10 hours
Kmode<-fit$K
fit$FeatureSal_out[[2]]
mcmc_samples<- fit$modeloutput_perchain$parameters.ECR.mcmc

sum_stats<-samplesSummary(mcmc_samples)


mu_rows <- grep("^mu\\[", rownames(sum_stats))
mu_summary <- sum_stats[mu_rows, ]

mu_matrix <- matrix(NA, nrow = Kmode, ncol = 18)

for(k in 1:Kmode) {
  for(j in 1:18) {
    row_name <- paste0("mu[", k, ", ", j, "]")
    if(row_name %in% rownames(mu_summary)) {
      mu_matrix[k, j] <- mu_summary[row_name, "Mean"]
    }
  }
}
