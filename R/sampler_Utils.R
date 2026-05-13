#Author: Carmen Rodriguez
#Created: 7/15/2025
#Last revised: 12/14/25
###REQUIRED PACKAGES:nimble,label.switching 

#===============================NIMBLE syntax===================================================
#Elements i need for nimble implementation
#1. Multivariate beta pdf  reparameterized using the mean and concentration parameter
#2. Complete data likelihood
#3. Priors for all my parameters; essentially priors for mu and lambda

#'@description   Multivariate Beta Distribution PDF
#https://r-nimble.org/cheatsheets/NimbleCheatSheet.pdf
#The multivariate beta distribution is not defined in NIMBLE;therefore, we will use nimbleFunction() to define a new distribution for use in NIMBLE
#See chapter 12 of the NIMBLE Manual
#A few notes:
#*double(0) argument and returnType(double(0)) establish that the input and output will both be zero-dimensional (scalar) numbers.
#*double(1): 1-dimension
#*To provide a user-defined distribution, you need to define density ("d") and optionally simulation (r)
#*NOTE:You need to provide the simulation (‘r’) function if any algorithm used with a model that uses the distribution needs to simulate from the distribution

#'@param x p-dimensional vector, (xi1, xi2,...,xip)--> row of X
#'@param mu  p-dimensional vector containing the expected value of feature xj
#'@param lambda concentration parameter remains constant across all dimensions, scalar
#'@param log for log of the PDF; default is not log
#'@return Value of the Multivariate Beta PDF for a single row 

dmultibeta <- nimbleFunction(
  run = function(x = double(1), mu = double(1), lambda = double(0), log = logical(0, default = 0)) {
    returnType(double())
    p <- length(x)
    delta <- mu / (1 - mu)
    sum_delta <- sum(delta)
    
    #Constant
    log_num <- lgamma(lambda * (1 + sum_delta))
    log_denom <- lgamma(lambda) + sum(lgamma(lambda * delta))
    
    # Numerator of the likelihood: ∏ x_j^(κ δ_j - 1)
    log_prod_x <- sum((delta * lambda - 1) * log(x) - (delta * lambda + 1) * log(1 - x))
    
    # Denominator: (1 + sum x / (1 - x))^(κ * (1 + sum δ))
    log_denom_x <- (lambda * (1 + sum_delta)) * log(1 + sum(x / (1 - x)))
    
    log_density <- log_num - log_denom + log_prod_x - log_denom_x
    
    if (log) return(log_density) else return(exp(log_density))
  }
)

#'@description   Multivariate Beta Distribution Random Generation
#'In this function, i follow the set up from Olkin & Liu
#'@param n number of observations. If length(n) > 1, the length is taken to be the number required
#'@param mu  p-dimensional vector containing the expected value of feature xj
#'@param lambda concentration parameter remains constant across all dimensions
#'@return  vector x ~ MBeta(mu,lambda)
#'
rmultibeta <- nimbleFunction(
  run = function(n = integer(0), mu = double(1), lambda = double(0)) {
    returnType(double(1))
    p <- length(mu)
    if (n != 1) print("Warning: only n = 1 is supported")
    
    delta <- mu / (1 - mu)
    alpha <- lambda  * delta
    #print(alpha)
    beta <- lambda 
    #print(beta)
    
    y <- numeric(p, init = TRUE)
    for (j in 1:p) {
      y[j] <- rgamma(1, shape = alpha[j], rate = 1)
    }
    u <- rgamma(1, shape = beta, rate = 1)
    
    x <- y / (y + u)  
    return(x)
  }
)


# dlogsumexp <- nimbleFunction(
#   run = function(loglik = double(1),
#                  log = integer(0, default = 0)) {
#     returnType(double(0))
#     
#     # numerically stable log-sum-exp
#     m <- max(loglik)
#     val <- m + log(sum(exp(loglik - m)))
#     
#     if (log) return(val) else return(exp(val))
#   }
# )


#nimbleCode ( add here because it's kinda long)--NEED TO ADD DESCRIPTION!
mbmmSF_buildnimbleCode<-function(Kmax,g, h, a_mu0, b_mu0 ,a_mu1, b_mu1,
                                 a_rho, b_rho, e0, priorOnAlpha){
  if (priorOnAlpha == "no"){
    modelcodeFS <- nimbleCode({
      #----Priors for model parameters----#
      
      # Mixture weights
      alpha[1:Kmax] <- rep(e0, Kmax) 
      pi[1:Kmax] ~ ddirch(alpha = alpha[1:Kmax])
      
      #Concentration parameters
      for (k in 1:Kmax){lambda[k] ~ dgamma(shape = g, scale=h)}
      
      #--Feature Saliency---# 
      #Feature Means 
      #Probability  for inclusion
      rho ~ dbeta(shape1 = a_rho, shape2 = b_rho)
      for (j in 1:p){
        #Variable relevance indicator
        s[j] ~ dbern(prob = rho) 
        # Binary Gibbs sampler: sampler_binary()
        
        #Mu for irrelevant variables
        mu0[j] ~ dbeta(shape1 = a_mu0, shape2 = b_mu0)
        
        # Cluster-specific for relevant means
        for (k in 1:Kmax) { mu1[k,j] ~ dbeta(shape1 = a_mu1,shape2 =  b_mu1)}
      }
      
      #Define component-specific means using saliency indicator
      for (k in 1:Kmax) {
        for (j in 1:p){mu[k, j] <- (1 - s[j]) * mu0[j] + s[j] * mu1[k, j]}}
      
      #Track number of filled components, that way we can do a trace plot for this
      #for (k in 1:Kmax) {n_k[k] <- sum(z[1:n] == k)}
      
      # Compute number of non-empty components
      #K0 <- Kmax - sum(n_k[1:Kmax] == 0)
      
      #--Complete Data Likelihood--#
      for (i in 1:n){
        z[i] ~ dcat(pi[1:Kmax])
        x[i, 1:p] ~ dmultibeta(mu[z[i], 1:p], lambda[z[i]])}
    })
    
  }else if (priorOnAlpha == "gam_05_05"){
    modelcodeFS2 <- nimbleCode({
      #----Priors for model parameters----#
      
      # Mixture weights
      alpha ~ dgamma(0.5, 0.5)
      for (k in 1:Kmax) {alpha_vec[k] <- alpha / Kmax}
      pi[1:Kmax] ~ ddirch(alpha =alpha_vec[1:Kmax])
      
      #Concentration parameters
      for (k in 1:Kmax){lambda[k] ~ dgamma(shape = g, scale=h)}
      
      #--Feature Saliency---# 
      #Feature Means 
      rho ~ dbeta(shape1 = a_rho, shape2 = b_rho)
      for (j in 1:p){
        s[j] ~ dbern(prob = rho) 
        mu0[j] ~ dbeta(shape1 = a_mu0, shape2 = b_mu0)
        for (k in 1:Kmax){mu1[k,j] ~ dbeta(shape1 = a_mu1,shape2 =  b_mu1)}
      }
      for (k in 1:Kmax) {
        for (j in 1:p){mu[k, j] <- (1 - s[j]) * mu0[j] + s[j] * mu1[k, j]}}
      
      #for (k in 1:Kmax) {n_k[k] <- sum(z[1:n] == k)}
      #K0 <- Kmax - sum(n_k[1:Kmax] == 0)
      
      #--Complete Data Likelihood--#
      for (i in 1:n) {
        z[i] ~ dcat(pi[1:Kmax])
        x[i, 1:p] ~ dmultibeta(mu[z[i], 1:p], lambda[z[i]])}
    })
    
  }else if (priorOnAlpha == "gam_1_2"){
    modelcodeFS2 <- nimbleCode({
      #----Priors for model parameters----#
      
      # Mixture weights
      alpha ~ dgamma(1, 2)
      for (k in 1:Kmax){alpha_vec[k] <- alpha / Kmax}
      pi[1:Kmax] ~ ddirch(alpha =alpha_vec[1:Kmax])
      
      #Concentration parameters
      for (k in 1:Kmax){lambda[k] ~ dgamma(shape = g, scale=h)}
      
      #--Feature Saliency---# 
      #Feature Means 
      rho ~ dbeta(shape1 = a_rho, shape2 = b_rho)
      for (j in 1:p){
        #Variable relevance indicator
        s[j] ~ dbern(prob = rho)
        mu0[j] ~ dbeta(shape1 = a_mu0, shape2 = b_mu0)
        for (k in 1:Kmax){mu1[k,j] ~ dbeta(shape1 = a_mu1,shape2 =  b_mu1)}
      }
      for (k in 1:Kmax){
        for (j in 1:p){mu[k, j] <- (1 - s[j]) * mu0[j] + s[j] * mu1[k, j]}}
      
      #for (k in 1:Kmax){n_k[k] <- sum(z[1:n] == k)}
      #K0 <- Kmax - sum(n_k[1:Kmax] == 0)
      
      #--Complete Data Likelihood--#
      for (i in 1:n){
        z[i] ~ dcat(pi[1:Kmax])
        x[i, 1:p] ~ dmultibeta(mu[z[i], 1:p], lambda[z[i]])}
    })
  }
  
}

#check<-mbmmSF_buildnimbleCode(Kmax=2,g=1, h=1, a_mu0=1, b_mu0=1 ,a_mu1=2, b_mu1=2,
#                       a_rho=1, b_rho=10, e0=1, priorOnAlpha = "no")


#create arrays for classification probs and label switching --CAN TRANSPOSE! FIX THIS FUNCTION: mcmc.pars1<-aperm(mcmc.pars, c(1,3,2))
#takes one chain at a time
#s: samples
#n: number of obs-- can remove this
#p: number of features
#K: number of components
#mcmc.pars  mxKxJ
#mcmc.pars1 mxJxK
mcmc_arraysv1<-function(s, n, p, K){
  m<- dim(s)[1]
  J<-p + 2 # muj parameter + lambda + pi
  
  # #Array mxKxJ-- for label switching
  # mcmc.pars<- array(data=NA, dim=c(m,K,J))
  # #pi
  # mcmc.pars[,,1]<- s[,paste0("pi[", 1:K, "]")]
  # #lambda
  # mcmc.pars[,,2]<-s[,paste0("lambda[", 1:K, "]")]
  # #mu
  # for (j in 1:p){
  #   mcmc.pars[,,2+j]<- s[, paste0("mu[", 1:K, ", ", j, "]")]
  # }
  #  dimnames(mcmc.pars) <- list(NULL, NULL, param = c("pi", "lambda", paste0("mu", 1:p)))
  
  #Array mxJxK for posterior probability computation
  mcmc.pars<- array(data=NA, dim=c(m,J,K))
  for (k in 1:K){
    mcmc.pars[,1,k]<-s[,paste0("pi[", k, "]")]
    mcmc.pars[,2,k]<-s[,paste0("lambda[", k, "]")]
    for (j in 1:p){
      mcmc.pars[,2+j,k]<- s[, paste0("mu[", k, ", ", j, "]")]
    }
  }
  dimnames(mcmc.pars) <- list(NULL, param = c("pi", "lambda", paste0("mu", 1:p)), component =  paste0("comp", 1:K))
  
  # Array mxKxJ-- for label switching
  mcmc.pars1<-aperm(mcmc.pars, c(1,3,2))
  dimnames(mcmc.pars1) <- list(NULL, NULL, param = c("pi", "lambda", paste0("mu", 1:p)))
  
  return(list("pp"=mcmc.pars, "ls"=mcmc.pars1))
}

mcmc_arrays<-function(mcmc_keep,p, Kmax,z){
  m<- nrow(mcmc_keep)
  J<-p + 3 #p  mu_j parameter + lambda + pi + nk
  
  # Array (mxK(max)x J)-- for label switching + the pruning and reordering
  mcmc.pars<- array(data=NA, dim=c(m,Kmax,J))
  #pi
  mcmc.pars[,,1]<- mcmc_keep[,paste0("pi[", 1:Kmax, "]")]
  #n_k
  #mcmc.pars[,,2]<- mcmc_keep[,paste0("n_k[", 1:Kmax, "]")]
  #Re-calculating n_k from z samples
  nk_matrix <- t(apply(z, 1, function(row) {
    # factor() ensures that even empty clusters get a count of 0
    as.numeric(table(factor(row, levels = 1:Kmax)))
  }))
  mcmc.pars[, , 2] <- nk_matrix
  

  #lambda
  mcmc.pars[,,3]<-mcmc_keep[,paste0("lambda[", 1:Kmax, "]")]
  #mu_j
  for (j in 1:p){
    mcmc.pars[,,3+j]<- mcmc_keep[, paste0("mu[", 1:Kmax, ", ", j, "]")]} #mu[,j]}
  dimnames(mcmc.pars) <- list(NULL, NULL,parameter = c("pi","n_k", "lambda",paste0("mu_", 1:p)))
  
  return(mcmc.pars)
}


reorder_and_prune_array <- function(mcmc_array, K_mode, Kmax, z, p) {
  m <- dim(mcmc_array)[1]
  J <- dim(mcmc_array)[3]
  n <- ncol(z)
  
  pruned <- array(NA, dim = c(m, K_mode, J))
  perms  <- matrix(NA, nrow = m, ncol = Kmax)
  z_pruned <- matrix(NA, nrow = m, ncol = n)
  
  for (t in 1:m) {
    nk_t <- as.numeric(mcmc_array[t, , "n_k"])
    nonempty <- which(nk_t > 0)
    stopifnot(length(nonempty) == K_mode)
    
    # permutation: active (alive) components first (van Havre)
    perm <- c(nonempty, setdiff(seq_len(Kmax), nonempty))
    perms[t, ] <- perm
    
    # relabel assignments for iteration t
    z_new <- match(z[t, ], perm)
    # drop labels > K_hat (these are empty components)
    z_new[z_new > K_mode] <- NA
    #print(z_new)
    z_pruned[t, ] <- z_new
   
     # reorder & prune all J parameter groups
    for (j in 1:J){
      tmp <- mcmc_array[t, perm, j]
      pruned[t, , j] <- tmp[1:K_mode]
    }
  }
  
  dimnames(pruned) <- list(NULL, NULL,parameter = c("pi",
"n_k","lambda",paste0("mu_", 1:p)))
  
  
  return(list("mcmcpruned" = pruned, "zpruned"= z_pruned))
}


#MAKE THIS INTO A NIMBLE FUNCTION SO IT GETS COMPILED!
#NIMBLE doesn’t  sample or monitor the posterior probabilities 
#Pr(zi=k|x_i) as part of the MCMC because z[i] is sampled from a categorical distribution with 
#probabilities defined by the likelihood × prior.
#But you can calculate those probabilities manually after sampling, using the saved MCMC samples.

#'@description function to calculate posterior probabilities by chain (from model output)
#' @param x data 
#' @param mcmc.pars mxJxK array of simulated MCMC samples
#' @param K number of components
#' @return mxnxK dimensional array of allocation probabilities
calc_zprobs_compiled <- nimbleFunction(
  run = function(x = double(2), 
                 pi = double(2),      # iterations x K
                 mu = double(3),      # iterations x K x p
                 lambda = double(2)   # iterations x K
  ) {
    returnType(double(3)) 
    m <- dim(pi)[1]
    K <- dim(pi)[2]
    n <- dim(x)[1]
    p <- dim(x)[2]
    z_probs <- array(0, dim = c(m, n, K))
    
    for(t in 1:m) {
      for(i in 1:n) {
        row_totals <- 0
        current_weights <- numeric(K)
        for(k in 1:K) {
          # Calls your existing dmultibeta C++ code directly
          log_dens <- dmultibeta(x[i, 1:p], mu[t, k, 1:p], lambda[t, k], log = TRUE)
          current_weights[k] <- exp(log(pi[t, k]) + log_dens)
          row_totals <- row_totals + current_weights[k]
        }
        if(row_totals > 0) {
          for(k in 1:K) { z_probs[t, i, k] <- current_weights[k] / row_totals }
        }
      }
    }
    return(z_probs)
  }
)


# Compile it into the project
C_calc_zprobs <- compileNimble(calc_zprobs_compiled)








# classificationprob<-function(x,n, K, m, J, mcmc.pars){
#   probs<-array(data = NA, dim=c(m,n,K)) #dimensionality of this array has to be  iteration x obs x K
#   for (iter in 1:m){
#     for (i in 1:n){
#       kdist <- c()
#       for (k in 1:K){
#         kdist[k]<- mcmc.pars[iter,"pi",k]*dmultibeta(x= x[i,], mu = mcmc.pars[iter,4:J,k] ,
#                                                      lambda = mcmc.pars[iter,"lambda",k])}
#       skdist <- sum(kdist)
#       probs[iter,i,]<-kdist/skdist
#     }}
#   return(probs)
# }
# 
# #nimbleFunctionList: PAGE 198 nimble user manual
# dmultibeta_class <- nimbleFunctionVirtual(
#   run = function(x = double(1),
#                  mu = double(1),
#                  lambda = double(0),
#                  log = integer(0, default = 0)) {
#     returnType(double())
#   }
# )
# classificationprob_nf <- nimbleFunction(
#   setup = function(dmbFun = dmultibeta_class){
#     dmb <- dmbFun    # compiled dmultibeta nimbleFunction
#   },
#   
#   run = function(x = double(2),      # n × p matrix
#                  n = integer(0),
#                  K = integer(0),
#                  m = integer(0),
#                  J = integer(0),
#                  mcmc_pars = double(3)) {
#     
#     returnType(double(3))   # m × n × K array
#     
#     probs <- array(0, c(m, n, K))
#     log_kdist <- numeric(K)
#     mu_vec <- numeric(J-3)
#     
#     for (iter in 1:m) {
#       for (i in 1:n) {
#         ## Compute unnormalized log-probabilities for each cluster
#         for (k in 1:K) {
#           # extract mu vector: indices 4:J
#           for (j in 4:J) {
#             mu_vec[j-3] <- mcmc_pars[iter, j, k]
#           }
#           
#           # log(pi_k)
#           log_pi_k <- log(mcmc_pars[iter, 1, k])
#           
#           # log-likelihood via dmultibeta (log=1)
#           log_dmb <- dmb$run(x[i, ], mu_vec, mcmc_pars[iter, 2, k], 1)
#           log_kdist[k] <- log_pi_k + log_dmb
#         }
#         
#         ## log-sum-exp normalizer
#         max_log <- max(log_kdist)
#         sum_exp <- 0.0
#         for (k in 1:K) {
#           sum_exp <- sum_exp + exp(log_kdist[k] - max_log)
#         }
#         log_norm <- max_log + log(sum_exp)
#         
#         ## normalized posterior classification prob
#         for (k in 1:K) {
#           probs[iter, i, k] <- exp(log_kdist[k] - log_norm)
#         }
#       }
#     }
#     
#     return(probs)
#   }
# )
# 



#Label Switching Algorithms as post-processing 

#Label Switching Algorithms 
#Supports 3 methods
#'@params zpivot p×n-dimensional array of pivot allocation vectors (ECR)
#'@params z m×n integer array of the latent allocation vectors generated from an MCMC algorithm.
#'@params K the number of mixture components.
#'@params zprobs. m×n×K dimensional array of allocation probabilities of the n observations among the 
#K mixture components, for each iteration  of the MCMC algorithm. 
#'@params method  subset of  c("STEPHENS", "ECR", "ECR-ITERATIVE-1"); can add multiple
#'@params mcmc.pars  mxKxJ array of the MCMC parameters
#'@return  an n dimensional vector of best clustering of the the observations for each method and relabeled parameters by method
#'mxKxJ
addressLabelSwitching<-function(method,zpivot,z,K, zprobs,mcmc.pars){
  ls <- label.switching(method = method,zpivot = zpivot, z = z, K = K, p = zprobs)
  L<-length(ls$permutations) #mxK array
  reordered.mcmc.pars<-list()
  for(l in 1:L){
    temp<-permute.mcmc(mcmc.pars, permutations =ls$permutations[[l]])
    reordered.mcmc.pars[[l]] <-temp$output}
  
  names(reordered.mcmc.pars)<-names(ls$permutations)
  #print(names(ls$permutations))
  #print(identical(mcmc.pars, reordered.mcmc.pars[["STEPHENS"]]))
  #print(identical(mcmc.pars, reordered.mcmc.pars[["ECR"]]))
  #print(identical(mcmc.pars, reordered.mcmc.pars[["ECR-ITERATIVE-1"]]))
  
  return(list("reordered.mcmc.pars"= reordered.mcmc.pars, "clusterMembershipPerMethod"=ls$clusters,
              "permutations" = ls$permutations))
}


###IN PROGRESS!!!
#Returns (for all chains):
# reordered.mcmc.pars: The aligned MCMC samples.
# classificationProbs: Posterior probabilities after relabeling.
# clusterMembership: Point estimate (e.g., MAP) cluster assignment per observation.
# all reordered parameters by the different methods so that we can use for later comparison


#mcmc.list: raw samples output from NIMBLE
#
# multiChainLabelSwitching<-function(x, mcmc.list, K){
#   n_chains <- length(mcmc.list)
#   #print(n_chains)
#   n <- nrow(x)  # number of observations
#   p<-ncol(x)   # number of features
#   results <- vector("list", n_chains)
#   
#   for (ch in seq_len(n_chains)){
#     # Extract data
#     this_chain <- as.matrix(mcmc.list[[ch]])
#     # print(class(this_chain))
#     #print(dim(this_chain))
#     mcmc.pars.arrays <- mcmc_arrays(s=this_chain,n = n,p=p,K = K)
#     #print(head(mcmc.pars.arrays))
#     pp<-mcmc.pars.arrays[["pp"]] #mxJxK
#     ls<-mcmc.pars.arrays[["ls"]] #mxKxJ
#     z <- this_chain[,tail(colnames(this_chain),n)]                 # m x n matrix of allocations
#     m <- dim(pp)[1]
#     J <- dim(pp)[2]
#     
#     # Compute initial (pre-relabeling) classification probabilities
#     zprobs_raw <- classificationprob(x = x, n=n, K = K, m=m, J=J, mcmc.pars = pp)
#     #print(head(zprobs_raw))
#     #pivot when method = ECR
#     #The pivot is selected by choosing a high-posterior density point, 
#     #such as the complete or non-complete Maximum A Posteriori (MAP) estimate. 
#     #most probable allocation to be used as the pivit labeling i.e. the one that best represents the posterior mode of cluster allocations
#     #MODE: the most frequently occurring configuration
#     mapindex<-which.max(table(apply(z, 1, paste, collapse = "-")))
#     #print(mapindex)
#     zpivot<- z[mapindex,]
#     
#     ls_result<-addressLabelSwitching(method = c("STEPHENS", "ECR", "ECR-ITERATIVE-1"), zpivot = zpivot, z= z,K = K, zprobs = zprobs_raw, mcmc.pars = ls)
#     #returns reordered samples and cluster membership per method and permutations
#     
#     # Pick one method’s relabeled parameters to proceed with (STEPHENS AS DEFAULT)-- similar a bayesbinmix los hicieron los dos
#     relabeled_pars <- ls_result$reordered.mcmc.pars[["STEPHENS"]] #mxKxJ
#     relabeled_pars_pp<- aperm(relabeled_pars, c(1,3,2)) #mxJxK
#     
#     # Recompute classification probabilities based on relabeled samples
#     relabeled_probs <- classificationprob(x = x, n=n, K = K, m=m, J=J, mcmc.pars = relabeled_pars_pp)
#     
#     #relabel cluster assignments in each MCMC iteration to make the labeling consistent across iterations.
#     #reordering allocations 
#     allocationsECR <- allocationsKL <- allocationsECR.ITERATIVE1 <-z
#     #print(names(ls_result$permutations))
#     #print(ls_result$permutations$"ECR"[1,])
#     for (i in 1:m){
#       #print(class(ls_result$permutations$"ECR"[i,]))
#       myPerm <- order(as.vector(ls_result$permutations$"ECR"[i,])) #need to make sure this is a vector
#       allocationsECR[i,] <- myPerm[z[i,]]
#       myPerm <- order(ls_result$permutations$"STEPHENS"[i,])
#       allocationsKL[i,] <- myPerm[z[i,]]
#       myPerm <- order(ls_result$permutations$"ECR-ITERATIVE-1"[i,])
#       allocationsECR.ITERATIVE1[i,] <- myPerm[z[i,]]
#     }
#     
#     #Reshape m × K × J into m × (K*J) Matrix to original form
#     param_names <- dimnames(relabeled_pars)[[3]]
#     param_mat <- matrix(NA, nrow = m, ncol = J * K)
#     col_names <- character(J * K)
#     for (k in 1:K){
#       idx <- ((k - 1) * J + 1):(k * J)
#       param_mat[, idx] <- relabeled_pars[, k, ]
#       col_names[idx]<-paste0("comp", k, "_", param_names)
#     }
#     colnames(param_mat) <- col_names
#     
#     
#     results[[ch]] <- list(
#       parameters.stephen.mcmc = mcmc(param_mat),
#       allocations.stephen.mcmc = mcmc(allocationsKL),
#       classificationProbabilities.stephen = relabeled_probs,
#       reordered_ALL = ls_result$reordered.mcmc.pars,
#       clusterMembershipPerMethod = ls_result$clusterMembershipPerMethod
#     )
#     
#     
#   } 
#   names(results) <- paste0("chain", seq_len(n_chains))
#   return(results)
# }


# Translated Poisson function to be use in NIMBLE- so maybe write as reg and a nimbleFunction
#' Translated Poisson Distribution
#' Density (`dtransPois`) and random generation (`rtransPois`) functions
#' for the translated Poisson distribution, where:
#' \deqn{K - 1 \sim \mathrm{Poisson}(s)}
#' so that \eqn{K = 1 + X}, \eqn{E[K] = s + 1}, and \eqn{Var[K] = s}.
#'
#' The support of `K` is \eqn{\{1, 2, 3, \ldots\}}.
#'
#' @param x Numeric vector of quantiles (for `dtransPois`).
#' @param s Positive numeric value for the Poisson rate parameter.
#' @param log Logical; if TRUE, probabilities `p` are returned on the log scale.
#' @param n Number of samples to generate (for `rtransPois`).
#'
#' @details
#' This distribution shifts the standard Poisson by +1. It can be used
#' to model count data that starts at 1 instead of 0.
# 
# transPois<- function(x, s, log= FALSE){
#   dpois(x - 1, lambda = s, log = log)
# }
# 
# dtransPois<- nimbleFunction(
#   run = function(x = double(0), s = double(0), log = integer(0, default = 0)) {
#     returnType(double(0))
#     if (log) 
#       return(dpois(x - 1, lambda = s, log = TRUE )) 
#     else 
#       return(dpois(x - 1, lambda = s))
#   }   
# )
# 
# rtransPois <- nimbleFunction(
#   run = function(n = integer(0), s = double(0)) {
#     returnType(double(0))
#     if (n != 1) print("rtransPois only allows n = 1; returning one sample.")
#     return(1 + rpois(1, lambda = s))
#   }
# )
    
#Quick check:
# s <- 3
# x <- 0:20
# k <- x + 1
# pk <- dpois(x, s)
# sum(k * pk)
# # Expected ≈ s + 1 = 4


 
  
