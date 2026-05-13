#Author: Carmen Rodriguez
#Helpful functions for data generations for simulations 
#Last revised: 7/15/2025
###REQUIRED PACKAGES: mclust,
library(mclust)
library(gtools)

#===============================Data generation======================================

#'@description   Multivariate Beta Distribution Random Generation
#'In this function, i follow the set up from Olkin & Liu
#'@param n number of observations. 
#'@param mu  p-dimensional vector containing the expected value of feature xj
#'@param lambda concentration parameter remains constant across all dimensions
#'@return  vector x ~ MBeta(mu,lambda)
#Generate data from a multivariate beta mixture model
#'@param
#'@param
#'@param
#'@return x nxp
generateData<-function(n,p, K, mu_mat, lambda_vec,p_vec){
  #counts <- round(p_vec * n)  # nk
  # # Create labels
 # labels <- rep(1:K, times = counts)
  #z_vec<-sample(labels)# shuffle
  z_vec <- sample(1:K, size = n, replace = TRUE, prob = p_vec)
  x <- matrix(NA, nrow = n, ncol = p)
  for (i in 1:n) {
    k <- z_vec[i]
    x[i, ] <- rmultibeta_other(n= 1, mu = mu_mat[k, ], lambda = lambda_vec[k])
  }
 return(list(x,z_vec))}


# generateData<-function(n,p, K, mu_mat, lambda_vec,p_vec){
#   alloc_vec<-rcat(n = n, prob = p_vec)
#   x <- matrix(NA, nrow = n, ncol = p)
#   for (i in 1:n) {
#     k <- alloc_vec[i]
#     x[i, ] <- rmultibeta_other(n= 1, mu = mu_mat[k, ], lambda = lambda_vec[k])
#   }
#   return(list(x,alloc_vec))}


# RNG for vector x ~ MBeta(mu,lambda): generate random deviates
rmultibeta_other <- function(n = integer(0), mu = double(1), lambda = double(0)) {
   # return a matrix of draws: 
    p <- length(mu)
    # delta <- mu / (1 - mu)
    # alpha <- lambda * delta
    # beta <- lambda
    # x <- matrix(0, nrow = n, ncol = p)
    # for (i in 1:n) {
    #   y <- numeric()
    #   for (j in 1:p) {
    #     y[j] <- rgamma(1, shape = alpha[j], rate = 1)
    #   }
    #   u <- rgamma(1, shape = beta, rate = 1)
    #   x[i, 1:p] <- y / (y + u)
    # }
    # return(x)
    if (n != 1) print("Warning: only n = 1 is supported")
    
    delta <- mu / (1 - mu)
    alpha <- lambda  * delta
    #print(alpha)
    beta <- lambda 
    #print(beta)
    
    y <- numeric()
    for (j in 1:p) {
      y[j] <- rgamma(1, shape = alpha[j], rate = 1)
    }
    u <- rgamma(1, shape = beta, rate = 1)
    
    x <- y / (y + u)  
    return(x)
  }


#======================Helpers for simulations===============

get_mode <- function(x) {
  uniq_x <- unique(x)
  uniq_x[which.max(tabulate(match(x, uniq_x)))]
}

#scenariosdat = alloclist[[1]]: data with true params for all scenarios
#cid: case ID
get_true_params<- function(scenariosdat, dgp_id){
  row <- dplyr::filter(scenariosdat, dgp_id == !!dgp_id)
  
  true_pi <- if(dgp_id %in% paste0("case",1:4)) dplyr::select(row, pi1, pi2) else row$pivec[[1]]
  
  true_mu <- row$mu[[1]]
  
  true_lambda <- if(dgp_id %in% paste0("case",1:4)) c(row$lambda1, row$lambda2) else row$lambda_vec[[1]]
  
  return(list(true_pi = as.numeric(true_pi), true_mu = true_mu,
              true_lambda = true_lambda, numfeat = row$p, Ktrue = row$Ktrue)) 
}
#Returns a list with all the component specific parameters
#get_true_params(scenariosdat, "case1")

##Note: function adapted from Dr.Wu SWOLCA get_dist() function
# `get_dist` returns the element-wise distance between 'par1' and 'par2' 
# according to the distance metric specified in 'dist_type'
#'@description  Calculate the MSE and mean absolute bias

# Inputs:
#   par1: First object. Can be vector, matrix, array
#   par2: Second object. Can be vector, matrix, array
#   dist_type: String specifying the distance type. Default is "mean_abs". Can
# also be  "mean_sq"
# Output: Element-wise distance
get_dist <- function(par1, par2,dist_type = "mean_abs") {
  if (dist_type == "mean_abs") {  # Mean absolute error
    dist <- mean(abs(par1 - par2))
  }else if (dist_type == "mean_sq") {  # MSE
    dist <- mean((par1 - par2)^2)
  } else {
    stop("Error: dist_type must be 'mean_abs', or 'mean_sq' ")
  }
  return(dist)
}



#' @title Run one simulation replicate
#' @description Fits the preK_mvb_nimblev2 model to one dataset, compares estimated
#' and true parameters, and computes MSEs and ARI for allocations.
#' @param sid Name of the scenario ("case1", "case2", "case3" for now)
#' @param rid Index of the dataset replicate (R=1:100)
#' @param sim_data List of simulated datasets by scenario
#' @param true_alloc List of true allocations by scenario
#' @param scenariosdat Tibble or data frame of true parameter values
#' @param K Number of mixture components to fit
#' @param niter, nburnin, nchain, ... Additional arguments passed to preK_mvb_nimblev2
#' @param p.true  An optional vector of cluster mixing weights considered as the ground-truth useful for simulations to align 
#' the estimated component parameters
#' @return A list with data.frame with MSE, MAB and ARI for the replicate, and a list with the estimated parameter values


# run_1sim_preK<-function(sid, rid, sim_data, true_alloc, sceneriosdat,
#                         K = 2, niter = 10000, nburnin = 5000,
#                         a_mu=0.5, b_mu=0.5, a_lambda=5, b_lambda=2,e0=1, priorOnAlpha = "no"){
#   
#   
#   casedat<- sim_data[[sid]][[rid]]
#   alloc_true<- true_alloc[[sid]][[rid]]
#   true_params<- get_true_params(scenariosdat, sid)
#   p<-true_params$numfeat
#   
#   
#   #Fit model
#   fit <- preK_mvb_nimblev2(
#     x = casedat,
#     K = K,
#     niter = niter,
#     nburnin = nburnin,
#     nchain = 1,
#     LambdaPrior = "gamma",
#     a_mu = a_mu, b_mu = b_mu,
#     a_lambda = a_lambda, b_lambda = b_lambda,
#     e0 = e0, priorOnAlpha = priorOnAlpha,
#     p.true = true_params$true_pi
#   )
#   
#   
#   #Posterior summaries
#   sSum <- samplesSummary(fit$modeloutput_perchain$parameters.stephen.mcmc)
#   
#   
#   #Now extract each parameter: see simis_basecase.R as reference
#   pi_est<- as.numeric(sSum[paste0("pi[", 1:K, "]"), 1])
#   lambda_est <- as.numeric(sSum[paste0("lambda[", 1:K, "]"), 1])
#   mu_mat_est <- sapply(1:p, function(j) {
#     as.numeric(sSum[paste0("mu[", 1:K, ", ", j, "]"), 1])
#   })
#   alloc_est<-fit$modeloutput_perchain$clusterMembershipPerMethod[1,] #Allocation after Stephens method for label switching
#   
#   
#   #Compute MSE and Mean Absolute bias for each parameter
#   mse_pi<-get_dist(pi_est, true_params$true_pi, dist_type = "mean_sq")
#   mab_pi<-get_dist(pi_est, true_params$true_pi, dist_type = "mean_abs")
#   
#   mse_lambda<-get_dist(lambda_est, true_params$true_lambda, dist_type = "mean_sq")
#   mab_lambda<-get_dist(lambda_est, true_params$true_lambda, dist_type = "mean_abs")
#   
#   mse_mu<-get_dist(mu_mat_est, true_params$true_mu, dist_type = "mean_sq")
#   mab_mu<-get_dist(mu_mat_est, true_params$true_mu, dist_type = "mean_abs")
#   
#   
#   ARI_alloc<- mclust::adjustedRandIndex(alloc_est, alloc_true)
#   
#   output<- list(data.frame(
#     scenario = sid,
#     replicate = rid,
#     mse_pi = mse_pi,
#     mse_lambda = mse_lambda,
#     mse_mu = mse_mu,
#     mab_pi= mab_pi,
#     mab_lambda = mab_lambda,
#     mab_mu = mab_mu,
#     ARIval = ARI_alloc), 
#     est_params = list(pi_est, lambda_est, mu_mat_est,alloc_est))
#   
#   return(output)
# }



#' @title Run one simulation replicate
#' @description Fits the Ombmm_FeatSaliency model to one dataset, compares estimated
#' and true parameters, and computes MSEs and ARI for allocations.
#' @param cid Name of the scenario ("case1", "case2", "case3" "case4", "case5)
#' @param rid Index of the dataset replicate (R=1:100)
#' @param sim_data List of simulated datasets by scenario
#' @param true_alloc List of true allocations by scenario
#' @param scenariosdat Tibble or data frame of true parameter values
#' @param Kmax Number of mixture components to fit
#' @param niter, nburnin, nchain, ... Additional arguments passed to preK_mvb_nimblev2
#' @param p.true  An optional vector of cluster mixing weights considered as the ground-truth useful for simulations to align 
#' the estimated component parameters
#' @param priorOnAlpha
#' @return A list with data.frame with MSE, MAB and ARI for the replicate, and a list with the estimated parameter values


run_1sim_OmbmmSF<-function(cid, rid, sim_data, true_alloc, sceneriosdat,
                        Kmax = 10, niter = 10000, nburnin = 3000,
                        g =  5, h = 3, a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
                        a_rho = 1, b_rho = 20, e0=0.02, priorOnAlpha = "no"){
  
  
  casedat<- sim_data[[cid]][[rid]]
  alloc_true<- true_alloc[[cid]][[rid]]
  true_params<- get_true_params(scenariosdat, cid)
  p<-true_params$numfeat
  Ktrue<-true_params$Ktrue
  
  #Fit model
  fit <- Ombmm_FeatSaliency(
    x = casedat,
    K = Kmax,
    niter = niter,
    nburnin = nburnin,
    nchain = 1,
    init_list = NULL, 
    g =  g, h = h, a_mu0 = a_mu0, b_mu0 = b_mu0,a_mu1 = a_mu1, b_mu1 = b_mu1,
    a_rho = a_rho, b_rho = b_rho, e0=0.01, priorOnAlpha = priorOnAlpha,
    p.true = true_params$true_pi, plot_trace = FALSE, calcWAIC = FALSE
  )
  
  if (is.null(fit)) {
    return(list(
      scenario  = cid,
      replicate = rid,
      Kmode     = 1,
      status    = "one cluster"
    ))
  }
  
  
  #Posterior summaries
  sSum <- samplesSummary(fit$modeloutput_perchain$parameters.ECR.mcmc)
  Kmode<- fit$K

  
  #If Kmode != Ktrue, skip all MSE/MAB calculations
  # ---------------------------------------------
  if (Kmode != Ktrue) {
    return( output = list(
      scenario = cid,
      replicate = rid,
      Kmode = Kmode,
      Ktrue = Ktrue,
      PIP = as.vector(fit$FeatureSal_out[[2]]),
      samples1 = fit$samples_raw,
      nk = fit$nk_iter
    ))
  }else{
  #Now extract each parameter: see simis_basecase.R as reference
  pi_est<- as.numeric(sSum[paste0("pi[", 1:Kmode, "]"), 1])
  lambda_est <- as.numeric(sSum[paste0("lambda[", 1:Kmode, "]"), 1])
  mu_mat_est <- sapply(1:p, function(j) {
    as.numeric(sSum[paste0("mu[", 1:Kmode, ", ", j, "]"), 1])
  })
  alloc_est<-fit$modeloutput_perchain$clusterMembershipPerMethod[1,] #Allocation after ECR-ITERATIVE-1 method for label switching
  
  
  #Compute MSE and Mean Absolute bias for each parameter
  mse_pi<-get_dist(pi_est, true_params$true_pi, dist_type = "mean_sq")
  mab_pi<-get_dist(pi_est, true_params$true_pi, dist_type = "mean_abs")
  
  mse_lambda<-get_dist(lambda_est, true_params$true_lambda, dist_type = "mean_sq")
  mab_lambda<-get_dist(lambda_est, true_params$true_lambda, dist_type = "mean_abs")
  
  mse_mu<-get_dist(mu_mat_est, true_params$true_mu, dist_type = "mean_sq")
  mab_mu<-get_dist(mu_mat_est, true_params$true_mu, dist_type = "mean_abs")
  
  
  ARI_alloc<- mclust::adjustedRandIndex(alloc_est, alloc_true)
  
  
  output<- list(data.frame(
    scenario = cid,
    replicate = rid,
    mse_pi = mse_pi,
    mse_lambda = mse_lambda,
    mse_mu = mse_mu,
    mab_pi= mab_pi,
    mab_lambda = mab_lambda,
    mab_mu = mab_mu,
    ARIval = ARI_alloc), 
    est_params = list(pi_est, lambda_est, mu_mat_est,alloc_est),
    PIP = as.vector(fit$FeatureSal_out[[2]]),
    nk = fit$nk_iter)
  
  return(output)}
  
  
}


#===============================Initialization functions======================================

#'@description  Initialize parameters using K-means and Method of Moments estimation
#' @param x  matrix containing data nxp
#' @param K number of clusters
#' @return list containing:
#'           1. Vectorized shape parameters for each mixture component/cluster (length = D*M)
#'           2. mixing proportions
#'           3. cluster assignment from K-means
#' 
init_params_kmeans<- function(x, K, eps=1e-6){
  # Initialize values using k-means and method of moments
  km.init<-kmeans(x, centers = K)
  mu_mat<-km.init$centers #Kxp
  # ensure mu_mat is strictly between 0 and 1
  mu_mat <- pmin(pmax(mu_mat, eps), 1 - eps)
  
  #  Initial mixing proportions assignments 
  nk<- sapply(1:K, function(i) length(which(km.init$cluster == i))) #number of obs in each cluster based on k-means
  p_vec<- nk/sum(nk) # mixing proportion: # cluster/ total
  
  # Cluster assignment
  z_vec<- km.init$cluster
  
  return(list(mu_mat, p_vec, z_vec))
}

#Note that GMM via mclust, chooses the best model based on BIC
#https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html#clustering
init_params_gmm<- function(x){
  # Initialize values using k-means and method of moments
  mod<-mclust::Mclust(x)
 #means
  mu_mat<-t(mod$parameters$mean)
  # mixing proportions assignments 
  p_vec<- mod$parameters$pro
  # cluster assignment
  z_vec<-mod$classification
  return(list(mu_mat, p_vec, z_vec))
}


#2x2
#See plots document
# init_params_base_2by2<-function(n,datatype = c("a", "b", "c"), a, b, c){
#   if (datatype=="a"){
#     #overlapping components; 2x2
#     p_vec<-c(0.35,0.65)
#     mu_mat<- matrix(c() nrow = 2, ncol= 2,byrow = T)
#     
#     lambda_vec<-c(30,35)
#     z_vec<-sample(1:2, n, replace = TRUE, prob = p_vec)
#     
  # }else if(datatype == "b"){
  #   #moderately separated components; 2x2 (iii from Trianasari et al. 2021)
  #   p_vec<-c(0.45,0.55)
  #   mu_mat<- matrix(c(35/80,25/80,15/85,35/85), nrow = 2, ncol= 2,byrow = T)    
  #   lambda_vec<-c(20,35)
  #   z_vec<-sample(1:2, n, replace = TRUE, prob = p_vec)}
  # else{
  #   # well separated components
  #   mu_mat <- matrix(c(
  #     0.15, 0.25,   # Component 1 (left-skewed)
  #     0.85, 0.75    # Component 2 (right-skewed)
  #   ), nrow = 2, byrow = TRUE)
  #   p_vec <- c(0.5, 0.5)
  #   z_vec <- sample(1:2, size = n, replace = TRUE, prob = p_vec)
  #   lambda_vec<-c(8,8)
  #   
#   # }
# return(list(mu_mat, p_vec, z_vec,lambda_vec))
# }

#===============================Other important functions======================================

#' #'Multivariate Beta PDF- Single observation
#' #'@param xi p-dimensional vector, (xi1, xi2,...,xip)--> row of X
#' #'@param mu_j  p-dimensional vector containing the expected value of feature xij
#' #'@param lambda concentration parameter
#' #'@return Value of the Multivariate Beta PDF for a single row 
#' 
#' mvb.pdf.i<- function(xi, mu_j, lambda){
#'   p<-length(mu_j)
#'   alpha_j = (mu_j/ (1-mu_j))*lambda
#'   beta = lambda
#'   
#'   alpha_norm<- sum(alpha_j) #sum_j=1^p (delta_j) -- sum across ROWS 
#'   #the Constant
#'   numerator<- gamma(lambda + alpha_norm)  
#'   denom<- gamma(lambda)*prod(gamma(alpha_j))
#'   constant <- numerator/denom
#'   #Other components
#'   pn<- prod(xi^(alpha_j - 1))#-- this does element wise vectorization
#'   pd<-prod((1-xi)^(alpha_j + 1)) #product on denominator
#'   s<-(1 +  sum((xi)/(1-xi)))^{-(lambda + alpha_norm)} #sum
#'   return(constant*(pn/pd)*s)
#' }

#' #'Multivariate Beta Mixture model PDF for each component k
#' #'@param xi p-dimensional vector, (xi1, xi2,...,xip)--> row of X
#' #' @param mu Kxp matrix 
#' #' @param lambdabold K-dimensional vector
#' #'@param pi vector of mixing proportions/probabilities 
#' #'@return Value of the Multivariate Beta Mixture PDF for each component j, for data point Xi
#' mvbmm.pdf.k<-function(xi, mu, lambdabold, pi){
#'   K<-length(pi)
#'   sapply(1:K, function(k) p[k]*mvb.pdf.i(xi,mu[k,], lambdabold[k]))
#' }
#' 
#' #'Multivariate Beta Mixture model PDF - sum of all components
#' #'@param Xi D-dimensional vector, (xi1, xi2,...,xid)--> row of X
#' #'@param  matrix containing shape parameters for each mixture component/cluster (alpha_j)
#' #'@param pj vector of mixing proportions/probabilities 
#' #'@return Value of the Multivariate Beta Mixture PDF for a data point xi
#' 
#' mvbmm.pdf<-function(xi, mu, lambdabold, pi){
#'   sum(mvbmm.pdf.k(xi, mu, lambdabold, pi))
#' }
#' 
#' 
#' #'Complete data Log- likelihood of the Multivariate Beta Mixture Model 
#' #'@param X Input data. nxp matrix
#' #'@param Z  vector of length n with values in 1:K indicating component assignments
#' #'@param mu Kxp; each row mu[k,] is the mean vector for component k
#' #'@param lambdabold vector of length K; concentration parameter for each component
#' #'@param pi vector of mixing proportions/probabilities 
#' #'@param K  number of components
#' #'@return total likelihood
#' 
#' complete_lik_mvb<- function(x, z, mu, lambda, pi, K) {
#'   n <- nrow(x)
#'   p <- ncol(x)
#'   total_lik <- 1  # Initialize
#'   for (k in 1:K) {
#'     idx <- which(z == k)
#'     nk <- length(idx)
#'     if (nk == 0) next
#'     
#'     xk <- x[idx, , drop = FALSE] #observations in cluster k
#'     mu_k <- mu[k, ]                      # Row k of mu (mean vector)
#'     delta <- mu_k / (1 - mu_k)           # p dimensional
#'     sum_delta <- sum(delta)  #sum_(j=1)^p delta_kj
#'     
#'     # First term: pi_k^{n_k}
#'     comp_lik <- pi[k]^nk
#'     
#'     # Constant term across observations assigned to component k
#'     gamma_num <- gamma(lambda[k] * (1 + sum_delta))
#'     gamma_denom <- gamma(lambda[k]) * prod(gamma(lambda[k] * delta))
#'     const_term <- gamma_num / gamma_denom
#'     
#'     obs_prod <- 1
#'     for (i in 1:nk) {
#'       xi <- xk[i, ]
#'       pn<- prod(xi^(lambda[k] * delta - 1))
#'       pd<-prod((1-xi)^(lambda[k] * delta + 1)) 
#'       odds <- xi / (1 - xi)
#'       s<-(1 + sum(odds))^(-lambda[k] * (1 + sum_delta)) 
#'       
#'       obs_prod <- obs_prod * ((pn/pd)*s)
#'     }
#'     
#'     # Component likelihood contribution
#'     comp_lik <- comp_lik * const_term^nk * obs_prod
#'     total_lik <- total_lik * comp_lik
#'   }
#'   
#'   return(total_lik)
#' }

