## Examining the behavior of the MBeta PDF and the mixture model ##
library(ggplot2)
library(nimble)
library(GGally)
library(dplyr)
library(purrr)
library(tibble)
source("R/helpers.R")

#---Univariate reparameterized-----
set.seed(123)

# Function using mean (mu) and concentration param (lambda)
#mu:controls the center of the distribution (where the mass is located)
#lambda:controls the tightness of the distribution (how concentrated or diffuse it is around mu)
reparam_beta_df <- function(n, mu, lambda) {
  alpha <- mu*lambda
  beta <- (1 - mu) * lambda
  tibble(
    x = seq(0.001, 0.999, length.out = n),
    density = dbeta(x, alpha, beta),
    mu = mu,
    lambda = lambda
  )
}

# Fix mu, vary lambda;
n<-2000
mu_fixed <- 0.6
lambda_vals <- c( 2, 5, 10, 20, 35, 50)
#?map_dfr(): binds dataframes rowwise

reparam_df <- map_dfr(lambda_vals, ~reparam_beta_df(n,mu_fixed, .x))

ggplot(reparam_df, aes(x = x, y = density, color = factor(lambda))) +
  geom_line(linewidth = 1) +
  labs(title = expression(paste("Beta PDF with ", mu, "=0.6, varying ", lambda)),
       color = "lambda", y = "Density") +
  theme_minimal()


# Fix sample size, vary mean
lambda_fixed <- 35
means <- seq(0.1, 0.9, length.out = 5)
mean_df <- map_dfr(means, ~reparam_beta_df(n,.x, lambda_fixed))
ggplot(mean_df, aes(x = x, y = density, color = factor(mu))) +
  geom_line(size = 1) +
  labs(title = expression(paste("Beta PDF with ", lambda, "=35, varying ", mu)),
       color = "μ (mean)", y = "Density") +
  theme_minimal()



#----Visualize how the common W (via lambda) changes dependence/shape---- 
set.seed(2008)

#Bi-variate
n<-5000
mu  <- c(.2, .8)
lambdas <- c(2, 5, 10, 20, 35, 50)


rmultibeta_n <- function(n, mu, lambda) {
  # return a matrix of draws:
  p <- length(mu)
   delta <- mu / (1 - mu)
  alpha <- lambda * delta
   beta <- lambda
  x <- matrix(0, nrow = n, ncol = p)
   for (i in 1:n) {
     y <- numeric()
     for (j in 1:p) {
       y[j] <- rgamma(1, shape = alpha[j], rate = 1)
    }
   u <- rgamma(1, shape = beta, rate = 1)
   x[i, 1:p] <- y / (y + u)
   }
   return(x)
}


out <- lapply(lambdas, function(l){
      X <- rmultibeta_n(n = n, mu = mu, lambda = l)
      df <- as.data.frame(X)
      df$lambda <- factor(l)
      df
    })
D <- do.call(rbind, out)
  
# Marginals
ggplot(D, aes(x = V1)) + 
    geom_density() + 
  facet_wrap(~lambda, nrow = 1, scales = "free_y") +
    labs(title = "", x = "Variable")

ggplot(D, aes(x = V2)) + 
  geom_density() + 
  facet_wrap(~lambda, nrow = 1, scales = "free_y") +
  labs(title = "V2 marginal vs lambda", x = "V2")

# ggplot(D, aes(x = V3)) + 
#   geom_density() + 
#   facet_wrap(~lambda, nrow = 1, scales = "free_y") +
#   labs(title = "V3 marginal vs lambda", x = "V3")
#   

# Pairwise scatter 
p1<-ggplot(D, aes(x = V1, y = V2)) + 
      geom_point() + facet_wrap(~lambda, nrow = 1) +
      labs(title = "")

p2<-ggpairs(D, columns = c("V1","V2"), ggplot2::aes(colour = lambda))

  

#Increase number of features p =4
  
n <- 5000
mu1 <- c(0.37,0.87,0.084,0.6)
lambdas1 <- c(5, 10, 35, 50)
out1 <- lapply(lambdas1, function(l){
    X <- rmultibeta_n(n = n, mu = mu1, lambda = l)
    df <- as.data.frame(X)
    df$lambda <- factor(l)
    df
  })
D1 <- do.call(rbind, out1)
  
# Marginals
ggplot(D1, aes(x = V1)) + 
  geom_density() + 
  facet_wrap(~lambda, nrow = 1, scales = "free_y") +
  labs(title = "V1 marginal vs lambda", x = "V1")

ggplot(D1, aes(x = V2)) + 
  geom_density() + 
  facet_wrap(~lambda, nrow = 1, scales = "free_y") +
  labs(title = "V2 marginal vs lambda", x = "V2")

ggplot(D1, aes(x = V3)) +
 geom_density() +
facet_wrap(~lambda, nrow = 1, scales = "free_y") +
labs(title = "V3 marginal vs lambda", x = "V3")

#
ggplot(D1, aes(x = V4)) +
  geom_density() +
 facet_wrap(~lambda, nrow = 1, scales = "free_y") +
 labs(title = "V4 marginal vs lambda", x = "V4")


 
p3<-ggpairs(D1, columns = c("V1","V2", "V3", "V4"), ggplot2::aes(colour = lambda))

  
  
## Mixture model
set.seed(123)
n <- 2500
p <- 4
K <- 2

# cluster parameters
mu_mat <- matrix(c(0.25, 0.65, 0.15, 0.4,
  0.75, 0.35, 0.55, 0.7), nrow = K, byrow = TRUE)

lambda_k1 <- c(15, 35)  # low vs high concentration
lambda_k2 <- c(10, 5) 
pi <- c(0.5, 0.5)         # equal mixing weights

# sample cluster memberships
z <- sample(1:K, n, replace = TRUE, prob = pi)

# generate data
x1 <- matrix(NA, n, p)
for (i in 1:n) {
  k <- z[i]
  x1[i, ] <- rmultibeta_other(n = 1, mu = mu_mat[k,], lambda = lambda_k1[k])
}


x2 <- matrix(NA, n, p)
for (i in 1:n) {
  k <- z[i]
  x2[i, ] <- rmultibeta_other(n = 1, mu = mu_mat[k,], lambda = lambda_k2[k])
}


Dmix <- as.data.frame(x1)
Dmix1<-as.data.frame(x2)
Dmix$Cluster <- factor(z)
Dmix1$Cluster <- factor(z)

# visualize pairwise relationships
p4<-ggpairs(Dmix, columns = 1:p,
        ggplot2::aes(color = Cluster, alpha = 0.5),
        upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 13)

p5<-ggpairs(Dmix1, columns = 1:p,
        ggplot2::aes(color = Cluster, alpha = 0.5),
        upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.3))) +
  theme_bw(base_size = 13)



pdf("simulations/figures/MbetaMM_behavior_plots.pdf")
p1
p2
p3
p4
p5
dev.off()

png("simulations/figures/Mbeta_lambda_behavior.png", width = 900, height = 300)
p1
dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#----Test whether a 1-component mixture is correctly recovered----
# Define a few scenarios
scenarios <- list(
  easy      = list(mu = c(0.2, 0.5, 0.8), lambda = 15),   # tight, easy to recover
  moderate  = list(mu = c(0.25, 0.55, 0.75), lambda = 5), # typical
  diffuse   = list(mu = c(0.3, 0.5, 0.7), lambda = 1.5),  # broad, challenging
  boundary  = list(mu = c(0.05, 0.5, 0.9), lambda = 6),   # with near-edges
  higherdim = list(mu = c(0.15, 0.3, 0.5, 0.7, 0.85), lambda = 6) # p=5 case
)

simulate_MB <- function(n, mu, lambda) {
  out <- rmultibeta_n(n = n, mu = mu, lambda = lambda)
  colnames(out) <- paste0("x", seq_along(mu))
  as.data.frame(out)
}

set.seed(42)
data_easy <- simulate_MB(n = 600, mu = scenarios$easy$mu, lambda = scenarios$easy$lambda)
head(data_easy)

ggplot(data_easy, aes(x = x1, y = x2)) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Scatter of first two dimensions", x = "x1", y = "x2")

plot_ly(data_easy, x = ~x1, y = ~x2, z = ~x3)
 

dat_mod<-simulate_MB(n = 1000, mu = scenarios$moderate$mu, lambda = scenarios$moderate$lambda)
head(dat_mod)

plot_ly(dat_mod, x = ~x1, y = ~x2, z = ~x3)
ggplot(dat_mod, aes(x = x1, y = x2)) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Scatter of first two dimensions", x = "x1", y = "x2")



##1) Fit a model with K = 1
# Dimensions
x<-as.matrix(dat_mod)
n <- nrow(x)
p <- ncol(x)



inits <- list(
  z = rep(1,1000),
  pi = 1,
  lambda = 5,
  mu = matrix(runif(1 * 3, 0.2, 0.8), nrow = 1)
)


# Default prior hyper parameters
constants <- list(n = n, p = p, K = 1,
                  a_lambda =  2, 
                  b_lambda = 1,
                  a_mu = 2, 
                  b_mu = 2,
                  alpha_pi = c(1))


# nimbleCode
#REMEMBER: No runtime conditionals (if) inside nimbleCode — everything is static
modelcode <- nimbleCode({
  lambda ~ dgamma(shape = a_lambda, scale = b_lambda)
  for (j in 1:p) { mu[j] ~ dbeta(shape1 = a_mu, shape2 = b_mu) }
  for (i in 1:n) { x[i, 1:p] ~ dmultibeta(mu[1:p], lambda) }
})


# Build model
mvb_model <- nimbleModel(
  code = modelcode,
  data = list(x = x),
  inits = list(lambda = 5, mu = runif(p, 0.2, 0.8)),
  constants = list(n = n, p = p, a_lambda = 2, b_lambda = 1, a_mu = 2, b_mu = 2),
  dimensions = list(mu = p, x = c(n, p))
)


# Configure and build MCMC
conf <- configureMCMC(mvb_model, monitors = c( "lambda", "mu"))
mcmc <- buildMCMC(conf, print = TRUE)


#Compile the NIMBLE model and the built MCMC algorithm into C++ for efficient execution.
Cmodel<- compileNimble(mvb_model)
Cmcmc <- compileNimble(mcmc,project = mvb_model)


# Run MCMC
samples <- runMCMC(Cmcmc, niter = 1000, nburnin = 500, nchains = 1,
                   setSeed = TRUE, samplesAsCodaMCMC = TRUE)



samplesSummary(samples)


##2) Using K > 1 and a sparse Dirichlet
test<-preK_mvb_nimble(x=x, K = 2, niter = 1000, nburnin = 500, nchain= 1, 
                a_lambda=2, b_lambda=1, a_mu=2, b_mu=2, e0=1/20)

head(test$modeloutput_perchain$classificationProbabilities.stephen)


test1<-preK_mvb_nimble(x=x, K = 3, niter = 1000, nburnin = 500, nchain= 1, 
                       a_lambda=2, b_lambda=1, a_mu=2, b_mu=2, e0=1/20)

head(test1$modeloutput_perchain$classificationProbabilities.stephen)



