#SAMBA-Mix Model
#Bayesian Overfitted multivariate beta mixture model with feature saliency
library(nimble)
library(coda)
library(patchwork)
library(bayesplot)
library(label.switching)
#Auxiliary functions
source("~/bayesmbmm/R/sampler_Utils.R")

#'@param x A numeric matrix of the data to be clustered. Rows are observations, columns are features. 
#' Must be of storage mode 'double'.
#'@param Kmax Integer; the maximum number of components in the overfitted mixture model.
#'@param niter A numeric scalar; specifying the number of MCMC iterations to run.
#'@param nburnin A numeric scalar; specifying the number of burn in samples to discard from each chain
#'@param nchain A numeric scalar; specifying the number of MCMC chains to execute
#'@param init_list initial values for each chain ; this can also be a function to generate initial values for each chain. 
#'Several functions are provided in sampler_Utils.R.
#'@param g,h hyperparameters for the gamma prior on lambda -- gamma( shape, scale)
#'@param a_mu0,b_mu0 hyperparameters for the beta prior on mu0: ADD MORE EXPLANATION
#'@param a_mu1,b_mu1 hyperparameters for the beta prior on mu1: ADD MORE EXPLANATION
#'@param priorOnAlpha  Character string; specifies the hyperprior on the Dirichlet concentration. 
#' Options include "no" (fixed e0), "gam_05_05", or "gam_1_2".
#'@param e0 Fixed scalar for Dirichlet prior weights when priorOnAlpha is "no". 
#' Small values (e.g., 0.01) induce sparsity (empty components).
#'@param a_rho,b_rho, hyperparameters for the beta prior on rho (prior inclusion probability (shared across variables))
#'@param p.true Optional numeric vector of ground-truth mixing weights for simulation performance metrics.
#'@param plot_trace Logical; if TRUE and calcWAIC is FALSE, generates trace plots for 
#' K0, pi, and cluster statistics.
#'@param calcWAIC  Logical; if TRUE, the function calculates the Marginal WAIC by 
#' integrating out the latent indicators (z). Note: This disables posterior 
#' tracking of z to ensure numerical stability in C++. See paper by Merkle et al. 2019, titled: Bayesian Comparison of Latent Variable Models: Conditional Versus Marginal Likelihoods
#'@param param_monitor Character vector of node names to monitor. Defaults to 
#' global parameters. Do not include "z" if calcWAIC is TRUE.

#'@return  list containing: 
#'WAIC: (Only if calcWAIC=TRUE) A nimbleList containing the WAIC, lppd, and pWAIC penalty.
#Only if calcWAIC=FALSE):
#' 1. nimblemodel specification
#' 2. nimblemodel configuration
#' 3. Raw MCMC samples (coda::mcmc object)
#' 4. Raw MCMC samples where K0 == K_mode
#' 5. K_mode
#' 6.FeatureSal_out list: Contains samples for the feature saliency indicators and Posterior inclusion probabilities (PIP)((p(s[j] =1)))
#' 7. modeloutput_perchain  list: (ONLY ONE CHAIN)
#'parameters.ECR.mcmc: object of class mcmc (see coda package) containing the simulated values (after burn-in)
#'allocations.ECR.mcmc:object of class mcmc (see coda package) containing the simulated values (after burn-in) of allocation variables
#'classificationProbabilities.ECR:   estimated posterior probabilities of belonging to a cluster from the multivariate beta mixture model.
#'reordered_ALL: reordered the MCMC samples according to ECR, STEPHENS and ECR-ITERATIVE-1 relabeling algorithms
#'clusterMembershipPerMethod: data frame of the most probable allocation of each observation after 
#'reordering the MCMC sample which corresponds to the most probable number of clusters according to  ECR-ITERATIVE-1 methods.
#' 8. chainInfo: number of mcmc cycles, burn-in period, and number of chains (1)


#default for rho induces sparsity: # ~ 0.0476
#Gamma(0.05, 0.05)-- very sparse same with the default e0= 0.1
Ombmm_FeatSaliency<- function(x, Kmax = 10, niter = 10000, nburnin = NULL, nchain= 1,init_list = NULL, 
                              g =  5, h = 3, a_mu0 = 1, b_mu0 = 1,a_mu1 = 2, b_mu1 = 2,
                              a_rho = 1, b_rho = 20, e0=0.01, priorOnAlpha ="no",
                              p.true= NULL, plot_trace = FALSE, calcWAIC = FALSE,
                              param_monitor =  c("s", "pi", "mu0", "mu1", "mu", "rho", "lambda", "z")){
  
  
  priorOnAlpha <- match.arg(priorOnAlpha, choices = c("no", "gam_05_05", "gam_1_2"))
  

  # Dimensions
  n <- nrow(x)
  p <- ncol(x)
  
  
  
  # nimbleCode definition: 
  # Build nimbleCode for different cases if hyperprior for dirichlet present
  modelcodeFS<-mbmmSF_buildnimbleCode(Kmax,g, h, a_mu0, b_mu0 ,a_mu1, b_mu1,
                                      a_rho, b_rho, e0, priorOnAlpha)
  

  #  Constants 
  constants <- list(
    n = n,
    p = p,
    Kmax = Kmax,
    e0 = e0,
    g =  g, h = h,
    a_mu0 = a_mu0, b_mu0 = b_mu0,
    a_mu1 = a_mu1, b_mu1 = b_mu1,
    a_rho = a_rho, b_rho = b_rho)
  
  
  # Data list:
  data <- list(x = x)
  
  # Default initial values if not provided
  
  if (is.null(init_list)){
    pi_init <- rep(1/Kmax, Kmax)
    pi_init <- pi_init / sum(pi_init) #Normalize
    init_list <- list(
      z= sample(1:Kmax, size = n, replace = TRUE),
      s = rep(1,p), 
      pi=pi_init,
      lambda = rgamma(Kmax, 9, 1),
      mu0 = rep(0.5, p),
      mu1 = matrix(0.5, nrow = Kmax, ncol = p))
  }
  
  # Build model: Configure and build MCMC
  #Select samplers and add monitors
  modelFS <- nimbleModel(code = modelcodeFS  , constants = constants, data = data, inits = init_list)
  #CmodelFS<- compileNimble(modelFS)
  
  #print(modelFS$isData("z"))
  #print(modelFS$getNodeNames(stochOnly = TRUE, includeData = FALSE))
  
  CmodelFS<- compileNimble(modelFS)

  # Note: If WAIC is being computed enableWAIC must be TRUE here--see main function inputs
  # Setup Marginal WAIC if calcWAIC = TRUE--> see Merkle et al. 2019
  z_nodes <- modelFS$expandNodeNames("z")
  conf <- configureMCMC(modelFS, monitors = param_monitor, enableWAIC = calcWAIC,
                        controlWAIC = if(calcWAIC) list(
                          marginalizeNodes = z_nodes, 
                          niterMarginal = 100, thin = TRUE
                        ) else list()
  )
  
  
  Rmcmc <- buildMCMC(conf, print = TRUE)
  
  
  #Compile the NIMBLE model and the built MCMC algorithm into C++ for efficient execution.
  Cmcmc <- compileNimble(Rmcmc, project = modelFS, resetFunctions = TRUE) 
  
  # Run MCMC
  samples <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin, nchains = nchain,
                     setSeed = TRUE, samplesAsCodaMCMC = TRUE, WAIC = calcWAIC)
  
  

  if (calcWAIC) {
    # --- EARLY EXIT FOR WAIC RUN ---
    message("\nWAIC Run Complete.")
    
    # Ensure WAIC was actually calculated to avoid NULL errors
    if (!is.null(samples$WAIC)) {
      waic_val <- samples$WAIC$WAIC
      p_waic   <- samples$WAIC$pWAIC
      lppd_val <- samples$WAIC$lppd
      
      cat("Marginal WAIC:", waic_val, "\n")
      cat("pWAIC (Penalty):", p_waic, "\n")
      cat("lppd:", lppd_val, "\n")
      
      # Return only the WAIC object and exit the function
      return(samples$WAIC) 
    } else {
      stop("WAIC was requested but not found in the output samples.")
    }
  }
  
  # --- PROCEED ONLY IF calcWAIC = FALSE ---
   message("\nInference Run Complete. Processing posterior samples...")
   
  #Get K0 = Kmax - sum(n_k[1:Kmax] == 0) to get the posterior distribution of the number of occupied components

  #Extract z-samples from the MCMC output
  # NIMBLE names them z[1], z[2], ..., z[n]
  z_cols <- grep("^z\\[", colnames(samples))
  if (length(z_cols) == 0) {
    stop("No 'z' samples found. Ensure 'z' is in param_monitor for inference runs.")
  }
  
  z_samples <- as.matrix(samples[, z_cols])# dimensions: iterations × n

  # Calculate n_k for every iteration
  # This creates a matrix where each row is an iteration and each column is a cluster (1 to Kmax)
  nk_matrix <- t(apply(z_samples, 1, function(row) {
    table(factor(row, levels = 1:Kmax))
  }))
  colnames(nk_matrix) <- paste0("n_k[", 1:Kmax, "]")
  
  #Calculate K0 (number of non-empty clusters) for every iteration
  K0_vec <- apply(nk_matrix, 1, function(row) sum(row > 0))
  
  #  K0 probability data frame
  K0df <- data.frame(K0 = K0_vec)
  post_prob_df <- K0df %>%
    count(K0) %>%
    mutate(posterior_prob = n / sum(n))
  message("\nPosterior probabilities for K0:\n")
  print(post_prob_df)
  

  #Explore posterior probability of n_k,pi_k, K0 and s
  #https://cran.r-project.org/web/packages/bayesplot/vignettes/visual-mcmc-diagnostics.html
  #Print plots for the user to see (let this be controlled via an on/off argument)
  if (plot_trace){
    #Posterior Distribution of K0
    p1<-ggplot(post_prob_df, aes(x = factor(K0), y = posterior_prob)) +
      geom_col(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
      labs(
        x = "K0 (number of non-empty components)",
        y = "Posterior Probability",
        title = "Posterior Distribution of K0"
      ) + theme_minimal()
    
    #Trace plots- Raw samples
    #Mixing Weights. (pi) Trace
    pi_samples<-samples[,paste0("pi[", 1:Kmax, "]")]
    pi_df <- as.data.frame(pi_samples)
    pi_long <- pi_df %>% pivot_longer(everything(),names_to = "Component", values_to = "Prob")
    pi_long$Iteration <- rep(1:nrow(pi_df), times = Kmax)
    # Plot all traces in one panel
    p2<- ggplot(pi_long, aes(x = Iteration, y = Prob, color = Component)) +
      geom_line(alpha = 0.7) +
      labs(title = "Trace plots for mixing weights (pi)-raw; pi_Kmax", y = expression(pi[k]), x = "Iteration") +
      theme_minimal() + scale_y_reverse() 
    
    #https://adiradaniel.netlify.app/post/ggmultipane/
    nested <- (p1/p2) + plot_annotation(tag_levels = 'A') #add figure labels
    nested #view multi-panel figure
    print(nested) 
  }
  
  #Number of non-empty components: K_mode = mode(K0)
  K_mode<-get_mode(K0_vec)
  message("\nMost probable number of clusters:\n")
  print(K_mode)
  
  #IF MODEL FINDS KMODE = 1, STOP!
  if (K_mode < 2) {
    warning(" Kmode < 2. Returning NULL.")
    return(NULL)
  }
  
  #Iterations with non-empty components (i.e., in how many iterations is K0 == K_hat)
  idx_iter<- which(K0_vec == K_mode)
  #print(K0_vec)
  

  #Approach from van Havre et al., 2015: Reorder and prune conditional on K_hat
  #Pruning must only occur on iterations that already have exactly K_hat filled components.
  mcmc_keep <- as.matrix(samples[idx_iter, ])
  #mcmc.pars <- mcmc_arrays(mcmc_keep,p, Kmax) #mxK(max)x J
  z<- as.matrix(mcmc_keep[,paste0("z[", 1:n, "]")])  # m x n matrix of allocations
  #print(colnames(z))
  #print(dim(z))
  
  mcmc.pars <- mcmc_arrays(mcmc_keep,p, Kmax, z) #mxK(max)x J
  
  mcmc_vH <- reorder_and_prune_array(mcmc.pars, K_mode, Kmax,z,p)
  #print(head(mcmc_vH$zpruned))
  #print(head(mcmc_vH$mcmcpruned))
  mcmc.pruned<-mcmc_vH$mcmcpruned #m x K_mode x J
  z_pruned<- mcmc_vH$zpruned 
  
  m <- dim(mcmc.pruned)[1]
  J <- dim(mcmc.pruned)[3]
  
  # Compute initial (pre-relabeling) classification probabilities
  pp<-aperm(mcmc.pruned, c(1,3,2)) #swap so that mxJxK_mode
  #zprobs_raw <- classificationprob(x = x, n=n, K = K_mode, m=m, J=J, mcmc.pars = pp)
  
  x_matrix <- as.matrix(x)
  storage.mode(x_matrix) <- "double"
  
  zprobs_raw <- C_calc_zprobs(
    x = x_matrix, 
    pi = mcmc.pruned[, , 1], 
    lambda = mcmc.pruned[, , 3],
    mu = mcmc.pruned[, , 4:dim(mcmc.pruned)[3]]
  )
  
  message("\nLabel switching process...")
  
  # Label Switching via label.switching R package: ECR Iterative
  #pivot when method = ECR
  # Maximum A Posteriori (MAP) estimate. 
  mapindex<-which.max(table(apply(z_pruned, 1, paste, collapse = "-")))
  zpivot<- z_pruned[mapindex,]
  
  ls_result<-addressLabelSwitching(method = "ECR-ITERATIVE-1", zpivot = zpivot, z = z_pruned, K = K_mode, mcmc.pars = mcmc.pruned)
  # Pick one method’s relabeled parameters to proceed with ECR-Iterative-1
  relabeled_pars <- ls_result$reordered.mcmc.pars$`ECR-ITERATIVE-1` #mx K_mode x J
  #Dimensions of this are in order: pi, nk, lambda, mu_j (j=1,...,p)
  
  # For simulations where we know the truth to keep label order consistent for MSE calculation
  # Also, for using the coda package we convert the samples back to the m x (K_mode * J) matrix format
  
  if (is.null(p.true) == FALSE){ #for simulations where we know the truth
    all_perms <- permutations(n = K_mode, r = K_mode) 
    #print(all_perms)
    pars_final <- array(NA, dim = c(m, K_mode, J))
    for (t in 1:m) {
      pi_vec <- relabeled_pars[t, , 1]  # extract pi vector for iteration i
      
      # Find the permutation that minimizes squared distance to true pi
      best_perm_index <- which.min(apply(all_perms, 1, function(p) {
        sum((pi_vec[p] - p.true)^2)
      }))
      best_perm <- all_perms[best_perm_index, ]
      #print(best_perm)
      
      # Reorder the component labels for all parameter blocks
      pars_final[t, , ] <-relabeled_pars[t, best_perm, ]
    }
    #Reshape m × K × J into m × (K*J) Matrix to original form
    param_mat <- matrix(NA, nrow = m, ncol = J * K_mode)
    col_idx <- 1
    for (j in 1:J) {
      for (k in 1:K_mode) {
        param_mat[, col_idx] <- pars_final[, k, j]
        col_idx <- col_idx + 1
      }
    }
    # Reorder so that "pi[...]" columns come first
    #NOTE (for example): the order if K=2 &  p=3 is mu[1,1], mu[2,1], mu[1,2], mu[2,2], mu[1,3], mu[2,3]
    
    pi_names <- paste0("pi[", 1:K_mode, "]")
    n_k_names<- paste0("n_k[", 1:K_mode, "]")
    lambda_names<-paste0("lambda[", 1:K_mode, "]")
    mu_names<- as.vector(sapply(1:p, function(j){paste0("mu[", 1:K_mode, ", ", j, "]")}, simplify = TRUE))
    reordered_names <- c(pi_names,n_k_names, lambda_names,mu_names)
    colnames(param_mat) <- reordered_names
    
  }else{
    #Reshape m × K × J into m × (K*J) Matrix to original form
    param_mat <- matrix(NA, nrow = m, ncol = J * K_mode)
    col_idx <- 1
    for (j in 1:J) {
      for (k in 1:K_mode) {
        param_mat[, col_idx] <- relabeled_pars[, k, j]
        col_idx <- col_idx + 1
      }
    }
    # Reorder so that "pi[...]" columns come first
    pi_names <- paste0("pi[", 1:K_mode, "]")
    n_k_names<- paste0("n_k[", 1:K_mode, "]")
    lambda_names<-paste0("lambda[", 1:K_mode, "]")
    mu_names<- as.vector(sapply(1:p, function(j){paste0("mu[", 1:K_mode, ", ", j, "]")}, simplify = TRUE))
    reordered_names <- c(pi_names, n_k_names, lambda_names,mu_names)
    colnames(param_mat) <- reordered_names
    
  }
  
  
  #Classification probabilities--final
  #relabel cluster assignments in each MCMC iteration to make the labeling consistent across iterations.
  #reordering allocations 
  allocationsECR.ITERATIVE1 <-z_pruned
  
  #Similar to BayesBinMix R package
  MeanReorderedpMatrix <- array(data = 0, dim = c(n,K_mode))    # define object that will contain the classification probs
  for (t in 1:m){
    myPerm <- order(ls_result$permutations$"ECR-ITERATIVE-1"[t,])
    allocationsECR.ITERATIVE1[t,] <- myPerm[z_pruned[t,]]
    #Classification probability matrix
    ecrperm<- ls_result$permutations$"ECR-ITERATIVE-1"[t,]# this is the permutation of labels for iteration t according to ECR iterative algorithm
    MeanReorderedpMatrix <- MeanReorderedpMatrix + zprobs_raw[t, ,ecrperm]   # apply permutation to the columns of zprobs_raw for given iteration and add the permuted matrix to MeanReorderedpMatrix
  }
  
  MeanReorderedpMatrix <- MeanReorderedpMatrix/m  
  
  modelParamsSummaries<- list(parameters.ECR.mcmc =  mcmc(param_mat),
                              allocations.ECR.mcmc = mcmc( allocationsECR.ITERATIVE1),
                              classificationProbabilities.ECR = MeanReorderedpMatrix,
                              reordered_ALL = ls_result$reordered.mcmc.pars,
                              clusterMembershipPerMethod = ls_result$clusterMembershipPerMethod)
  
  #Feature Saliency Output
  s_samps<- mcmc_keep[, paste0("s[", 1:p, "]")]
  #Posterior inclusion probabilities (PIP)
  pip<-colMeans(s_samps)
  
  
  
  return(list(nimblemodel = modelFS, modelconfig = conf,
               samples_raw = samples, samples_keep = mcmc_keep,
               K = K_mode, nk_iter= nk_matrix, 
               FeatureSal_out = list(s_samps, pip),
               modeloutput_perchain =  modelParamsSummaries,
               chainInfo = c(niter, nburnin, nchain)))
  
}
