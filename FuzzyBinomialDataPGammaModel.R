# Coding up Bayes Logit 

#install.packages("emdbook")
#install.packages("mvnfast") # <- Updated in March-April 2020 to speed up sampling 
#install.packages("BayesLogit")
#install.packages("tidyr")
# install.packages("TruncatedDistributions", repos="http://R-Forge.R-project.org")

#require(emdbook)
library(BayesLogit)
#require(mvnfast)
require(tidyr)
require(Matrix)
require(truncnorm)
require(microbenchmark)
#require(TruncatedDistributions)

# Data_Frames: list of X - fixed effects matrix, A_1 - random effect matrix, ..., A_R - Rth random effect matrix 
# Y = number of sites length vector of vaccination counts per site
# N = number of sites length vector with total number of vaccine-eligible veterans at each site
# nsamps = number of posterior samples 
# warmup = burnin stage (discarded posterior samples, enough for the sampler to 'reach convergence')
# hyp = vector of prior rate and scale parameters for inverse gamma prior on random effect variances
# these two numbers are used for all hiearchical levels/random effects' variances
# mu_delt = number of site-length vector of hypothesized mean number of veterans having unaccounted for vaccinations
# sig_delt = number of site-length vector of hypothesized 'trust' one has in mu_delt values per site to reflect true delta
# smaller sig_delt implies higher trust in mu_delt 
# a_vec = lower bound options for delta draw, number of site-length vector of lower bounds 
# if NULL and 'underreport' = TRUE, then a_vec = 0 for all sites 
# underreport = TRUE indicates you are using the proposed method, underreport = FALSE indicates you are simply
# using hierarchical logistic regression in a bayesian framework, no imputation step 
# VA_Only = FALSE indicates the focus is on vaccination rates in general
# VA_Only = TRUE indicates the focus is on vaccination rates in the VA only. 

site_pgamma1 = function(Data_Frames,Y,N,nsamps=4000,warmup=500, hyp = c(0.01, 0.01),mu_delt = NULL, sig_delt = NULL,
                       a_vec= NULL, underreport = TRUE, VA_only = FALSE){
  # this is on the site level, so nrow(X) = number of sites 
  # fixed effects first in Data_Frames list
  a = hyp[1]; b = hyp[2]
  l_q <- length(Data_Frames)
  dim_q <- sapply(1:l_q, function(x) ncol(Data_Frames[[x]]))
  
  j = 1; storseq <- list()
  
  for(k in 1:l_q){
    dim1 <- dim_q[k]
    if(k == 1){seq1 <- j:dim1}else{
      seq1 <- j:sum(dim_q[1:k])
    }
    storseq[[k]] <- seq1
    j <- max(seq1)+1
    }
  
  tot_effs <- sum(dim_q)
  n_site <- W <- nrow(Data_Frames[[1]])
  if(is.null(a_vec)){a_vec = rep(0,n_site)}
  X <- do.call(cbind,Data_Frames)
  
  # random jitter in outcomes:
  if(underreport){
  if(VA_only){
    delta <- ceiling(rtruncnorm(W,mean=mu_delt,sd= sqrt(sig_delt),a = a_vec, b = (N-(Y+1))))
    N_Star <- N-delta # what we think if VA site only
    Y_Star = Y
  }else{
    delta <- ceiling(rtruncnorm(W,mean=mu_delt,sd= sqrt(sig_delt),a = a_vec, b = (N-(Y+1))))
    N_Star <- N
    Y_Star = Y + delta # what we think if we include all possible data sources
  }
  }else{
  Y_Star = Y; N_Star = N
  }
  kappa = (Y_Star-0.5*N_Star)
  # Storage object
  beta.samps <- lapply(1:l_q, function(y) matrix(NA,nrow=nsamps,ncol = dim_q[y]))
  if(l_q > 1){
  tau.samps <- lapply(1:(l_q-1), function(x) rep(NA,nsamps))
  tau = rep(0.5, l_q) 
  tau = tau[-1]
  }else{
    tau.samps <- NULL; 
    tau = NULL
  }
  
  #inits
  beta_vec <- rep(1, tot_effs)
  B_inv <- rep(0.01, dim_q[1])
  matlist <- unlist(c(B_inv, c(sapply(2:length(dim_q),function(x) rep(tau[x-1],dim_q[x])))))
  Big_Prec <- diag(matlist, nrow = length(matlist))
  
  # Gibbs Sampler  
  for(samp in 1:(nsamps+warmup)){
    # Sample alpha
    #print(paste0("Posterior run ", samp, 'started'))
    
    omega = rpg.devroye(n_site,N_Star,X%*%beta_vec)        # Update auxillary parameters
    
    V_Beta = chol2inv(chol(crossprod(X,(X*omega)) + Big_Prec))
    m_beta = as.matrix(V_Beta%*%crossprod(X,kappa))
                  # Mean vector
    beta_vec = (m_beta + t(chol(V_Beta))%*%rnorm(ncol(X)))# Unconstrained

    if(l_q > 1){
    #sample from posterior for tau
    tau = sapply(2:l_q, function(x){
      apart <- t(beta_vec[storseq[[x]]])%*%beta_vec[storseq[[x]]]/2;
      tauout <- rgamma(1, shape = a + length(storseq[[x]])/2, rate = apart + b);
    return(tauout)})
    tau <- c(NA,tau);
    matlist <- unlist(c(B_inv, c(sapply(2:length(dim_q),function(x) rep(tau[x],dim_q[x])))))
    Big_Prec <- diag(matlist, nrow = length(matlist))
    }

    # random jitter in outcomes:
    if(underreport){
    if(VA_only){
      delta <- ceiling(rtruncnorm(W,mean=mu_delt,sd= sqrt(sig_delt),a = a_vec, b = (N-Y-1)))
      N_Star <- N-delta # what we think if VA site only
      Y_Star = Y
    }else{
      delta <- ceiling(rtruncnorm(W,mean=mu_delt,sd= sqrt(sig_delt),a = a_vec, b = (N-Y-1)))
      N_Star <- N
      Y_Star = Y + delta # what we think if we include all possible data sources
    }
    # we need to restrict this fraction to be smaller than 1. 
    # this is how. 
    kappa = (Y_Star-0.5*N_Star)
    }
    # Save posterior sample 
    if(samp>warmup){
      for(k in 1:l_q){
        beta.samps[[k]][samp-warmup,] <- beta_vec[storseq[[k]]]
        if(k > 1){
          tau.samps[[k-1]][samp-warmup] <- tau[k]
        }}
    } 
  }
  return(list(Beta_Samps = beta.samps,Tau_samps=tau.samps)) 
}

# conversion of a 'nested' (I think) ICCs to variances in random effects
#input is a length 1 or length 2 vector of ICCs for two levels or three levels hierarchies
# computes 
convertICCtoVar <- function(icc){
  l_icc <- length(icc)
  if(l_icc == 1){
    sigma_site <- icc*(pi^2/3)/(1 - icc)
    sigma_out <- sigma_site
  }
  if(l_icc == 2){
    first = (icc[2] + icc[2]*icc[1])
    sigma_hcl = (pi^2/3)*first / (1-first)
    sigma_site = (((pi^2/3) + sigma_hcl)*icc[1])/(1-icc[1])
    sigma_out <- c(sigma_site, sigma_hcl)
  }
  return(sigma_out)
}

# this second model handles misclassification on the individual person level 
# assumes there is a prior probability of misclassification for each person 
# estimates it and imputes a 'true value' instead for the analysis. 
# but posterior probabilities of being mislabeled are based on original Y

# Data_Frames: list of X - fixed effects matrix, A_1 - random effect matrix, ..., A_R - Rth random effect matrix 
# Y = number of observations length vector of vaccination status
# nsamps = number of posterior samples 
# warmup = burnin stage (discarded posterior samples, enough for the sampler to 'reach convergence')
# hyp = vector of prior rate and scale parameters for inverse gamma prior on random effect variances
# these two numbers are used for all hiearchical levels/random effects' variances
# ab_T = two hyperparameter values for the prior probability that for an individual observation 
# we should leave this observation alone (T_i = 1) or estimate if its misclassified (T_i = 0)
# the current specifications are from Verdinelli and Wasserman (1991)
# ab_gam = two hyperparameter values for the prior probability of misclassification IN GENERAL in this data
# set. So think ab_gam[1]/(ab_gam[1] + ab_gam[2]) as the expected prior probability of misclassification
# right now this prior probability is specified to be on average 10% misclassification
# Question_IDs = the IDs of individuals we are unsure of their vaccination status
# should be a vector of numbers corresponding to row indices in data set of Data_Frames

RobustMisclassification = function(Data_Frames,Y,nsamps=4000,warmup=500, hyp = c(0.01, 0.01),
                                   ab_T = c(2,18), ab_gam = c(2,18), Question_IDs = NULL){
  # this is on the site level, so nrow(X) = number of sites 
  # fixed effects first in Data_Frames list
  a = hyp[1]; b = hyp[2] # hyper parameter specifications for variances
  a_T = ab_T[1]; b_T = ab_T[2] # prior specification on T - shared probability of correct classification
  a_gam = ab_gam[1]; b_gam = ab_gam[2] # prior specification for gamma
  LongN <- nrow(Data_Frames[[1]])
  # if we have known unsure veteran IDs
  if(!is.null(Question_IDs)){
  length_Q <- length(Question_IDs)
  }else{Question_IDs = 1:LongN}
  
  l_q <- length(Data_Frames)
  dim_q <- sapply(1:l_q, function(x) ncol(Data_Frames[[x]]))
  
  j = 1; storseq <- list()
  
  for(k in 1:l_q){
    dim1 <- dim_q[k]
    if(k == 1){seq1 <- j:dim1}else{
      seq1 <- j:sum(dim_q[1:k])
    }
    storseq[[k]] <- seq1
    j <- max(seq1)+1
  }
  
  tot_effs <- sum(dim_q)
  X <- do.call(cbind,Data_Frames)
  
  #initialize
  kappa <- Y - 1/2
  gammas <- .05
  
  # Storage object
  beta.samps <- lapply(1:l_q, function(y) matrix(NA,nrow=nsamps,ncol = dim_q[y]))
  if(l_q > 1){
    tau.samps <- lapply(1:(l_q-1), function(x) rep(NA,nsamps))
    tau = rep(0.5, l_q) 
    tau = tau[-1]
  }else{
    tau.samps <- NULL; 
    tau = NULL
  }
  delta.samps <- numeric(LongN)
  epsilons <- rep(.05, LongN)
  epsilon.samps <- numeric(length_Q)
  gamma.samps <- numeric(nsamps)
  if(!is.null(Question_IDs)){
    delta.samps <- numeric(length_Q)
    epsilon.samps <- numeric(length_Q)
    epsilons <- rep(.05, length_Q)
  }
  
  #inits
  beta_vec <- rep(1, tot_effs)
  B_inv <- rep(0.01, dim_q[1])
  matlist <- unlist(c(B_inv, c(sapply(2:length(dim_q),function(x) rep(tau[x-1],dim_q[x])))))
  Big_Prec <- diag(matlist, nrow = length(matlist))
  
  Y_old <- Y
  # Gibbs Sampler  
  for(samp in 1:(nsamps+warmup)){
    # Sample alpha
    print(paste0("Posterior run ", samp, ' started'))
    # individual level model 
    omega = rpg.devroye(LongN,1,X%*%beta_vec)        # Update auxillary parameters
    
    V_Beta = chol2inv(chol(crossprod(X,(X*omega)) + Big_Prec))
    m_beta = as.matrix(V_Beta%*%crossprod(X,kappa))
    # Mean vector
    beta_vec = (m_beta + t(chol(V_Beta))%*%rnorm(ncol(X)))# Unconstrained
    
    # expits_vec
    expits_vec <- c(expit(X[Question_IDs,]%*%beta_vec))
    # evaluate P_T: Ti is spike versus slab probability 
    P_T <- P_T_Eval(expits_vec, gammas, Y_old[Question_IDs], epsilons, a_T, b_T)
    # sample for T_i indicator 
    T_i <- rbinom(length_Q, 1, P_T)
    
    # sample switching probabilities
    epsilons <- P_Eps_Eval(expits_vec, a_T, b_T, Y_old[Question_IDs], length_Q,T_i)
    
    # sample for misclassification indicators delta
    deltas <- Eval_Deltas(expits_vec, epsilons, Y_old[Question_IDs], length_Q)
    # grab new projected Y_is
    Y_stars <- deltas[[2]]
    deltas <- deltas[[1]]
    
    # sample for shared probability of misclassification
    # beta distribution
    gammas <- Samp_Gamma(a_gam, b_gam, T_i, length_Q)
    
    Y[Question_IDs] <- Y_stars #working outcomes
    kappa <- Y - 1/2
    
    if(l_q > 1){
      #sample from posterior for tau
      tau = sapply(2:l_q, function(x){
        apart <- t(beta_vec[storseq[[x]]])%*%beta_vec[storseq[[x]]]/2;
        tauout <- rgamma(1, shape = a + length(storseq[[x]])/2, rate = apart + b);
        return(tauout)})
      tau <- c(NA,tau);
      matlist <- unlist(c(B_inv, c(sapply(2:length(dim_q),function(x) rep(tau[x],dim_q[x])))))
      Big_Prec <- diag(matlist, nrow = length(matlist))
    }
    
    # Save posterior sample 
    if(samp>warmup){
      for(k in 1:l_q){
        beta.samps[[k]][samp-warmup,] <- beta_vec[storseq[[k]]]
        if(k > 1){
          tau.samps[[k-1]][samp-warmup] <- tau[k]
        }}
      delta.samps <- deltas + delta.samps
      epsilon.samps <- epsilons + epsilon.samps
      gamma.samps[samp-warmup] <- gammas
    } 
  }
  return(list(Beta_Samps = beta.samps,Tau_samps=tau.samps, Delta_Samps = delta.samps/(nsamps-warmup), Gamma_Samps = gamma.samps,
              Epsilon_Samps = epsilon.samps/(nsamps-warmup))) 
}

expit <- function(x){exp(x)/(1+exp(x))}
# assistant functions for the individual based misclassification model
P_T_Eval <- function(expits_vec, gammas, Y, epsilons, c, d){
  numerator <- gammas*(expits_vec^Y)*((1-expits_vec)^(1-Y))
  
  ILFirst <- expits_vec*gamma(c + (1-Y))*gamma(d+Y)
  ILSecond <- (1-expits_vec)*gamma(c + Y)*gamma(d+ (1-Y))
  IntegrateL <- (ILFirst + ILSecond)/((c+d)*gamma(c)*gamma(d))
  
  second <- (1-gammas)*IntegrateL
  
  prob_Ti <- numerator/(numerator+second); 
  rm(numerator); rm(ILFirst); rm(ILSecond);rm(IntegrateL);rm(second)
  return(prob_Ti)
}

P_Eps_Eval <- function(expits_vec, c, d, Y, LongN, T_i){
  Z_i <- rbinom(LongN, 1, prob = expits_vec)
  eps_out <- Z_i*rbeta(LongN, c+(1-Y), d + Y) + ((1-Z_i) * rbeta(LongN, c+Y, d + (1-Y)))
  eps_out[T_i == 1] <- 0
  return(eps_out)
}

Eval_Deltas <- function(expits_vec, epsilons, Y, LongN){
  numerators <- epsilons*((1-expits_vec)^Y)*(expits_vec^(1-Y))
  denominator2 <- (1-epsilons)*(expits_vec^Y)*((1-expits_vec)^(1-Y))
  fractions <- numerators/(numerators + denominator2)
  deltas <- rbinom(LongN, 1, fractions)
  Y_star <- (1-deltas)*Y + (deltas)*(1-Y)
  return(list('Deltas' = deltas, 'Y_star' = Y_star))
}

Samp_Gamma <- function(c_gam, d_gam, T_i, LongN){
  first <- c_gam + sum(T_i)
  second <- d_gam + LongN - sum(T_i)
  outgam <- rbeta(1, first, second)
  return(outgam)
}

# simpler robust regression, shared error probability 
# distribution (misclassification epsilon)

RobustMisclassification_Simp = function(Data_Frames,Y,nsamps=4000,warmup=500, hyp = c(0.01, 0.01),
                                   ab_T = c(2,18), Question_IDs = NULL){
  # this is on the site level, so nrow(X) = number of sites 
  # fixed effects first in Data_Frames list
  a = hyp[1]; b = hyp[2] # hyper parameter specifications for variances
  a_T = ab_T[1]; b_T = ab_T[2] # prior specification on T - shared probability of correct classification
  LongN <- nrow(Data_Frames[[1]])
  # if we have known unsure veteran IDs
  if(!is.null(Question_IDs)){
    length_Q <- length(Question_IDs)
  }else{Question_IDs = 1:LongN}
  
  l_q <- length(Data_Frames)
  dim_q <- sapply(1:l_q, function(x) ncol(Data_Frames[[x]]))
  
  j = 1; storseq <- list()
  
  for(k in 1:l_q){
    dim1 <- dim_q[k]
    if(k == 1){seq1 <- j:dim1}else{
      seq1 <- j:sum(dim_q[1:k])
    }
    storseq[[k]] <- seq1
    j <- max(seq1)+1
  }
  
  tot_effs <- sum(dim_q)
  X <- do.call(cbind,Data_Frames)
  
  #initialize
  kappa <- Y - 1/2
  
  # Storage object
  beta.samps <- lapply(1:l_q, function(y) matrix(NA,nrow=nsamps,ncol = dim_q[y]))
  if(l_q > 1){
    tau.samps <- lapply(1:(l_q-1), function(x) rep(NA,nsamps))
    tau = rep(0.5, l_q) 
    tau = tau[-1]
  }else{
    tau.samps <- NULL; 
    tau = NULL
  }
  delta.samps <- numeric(LongN)
  deltas <- rep(0, LongN)
  epsilons <- .5
  epsilon.samps <- numeric(nsamps)
  if(!is.null(Question_IDs)){
    delta.samps <- numeric(length_Q)
    deltas <- rep(0, length_Q)
      }

  #inits
  beta_vec <- rep(1, tot_effs)
  B_inv <- rep(0.01, dim_q[1])
  matlist <- unlist(c(B_inv, c(sapply(2:length(dim_q),function(x) rep(tau[x-1],dim_q[x])))))
  Big_Prec <- diag(matlist, nrow = length(matlist))
  
  Y_old <- Y
  # Gibbs Sampler  
  for(samp in 1:(nsamps+warmup)){
    # Sample alpha
    print(paste0("Posterior run ", samp, ' started'))
    # individual level model 
    omega = rpg.devroye(LongN,1,X%*%beta_vec)        # Update auxillary parameters
    
    V_Beta = chol2inv(chol(crossprod(X,(X*omega)) + Big_Prec))
    m_beta = as.matrix(V_Beta%*%crossprod(X,kappa))
    # Mean vector
    beta_vec = (m_beta + t(chol(V_Beta))%*%rnorm(ncol(X)))# Unconstrained
    
    # expits_vec
    expits_vec <- c(expit(X[Question_IDs,]%*%beta_vec))
    
    # sample switching probabilities
    epsilons <- P_Eps_Eval_Simp(a_T, b_T, length_Q,deltas)
    
    # sample for misclassification indicators delta
    deltas <- Eval_Deltas(expits_vec, epsilons, Y_old[Question_IDs], length_Q)
    # grab new projected Y_is
    Y_stars <- deltas[[2]]
    deltas <- deltas[[1]]
    
    # sample for shared probability of misclassification
    # beta distribution
    #gammas <- Samp_Gamma(a_gam, b_gam, T_i, length_Q)
    
    Y[Question_IDs] <- Y_stars #working outcomes
    kappa <- Y - 1/2
    
    if(l_q > 1){
      #sample from posterior for tau
      tau = sapply(2:l_q, function(x){
        apart <- t(beta_vec[storseq[[x]]])%*%beta_vec[storseq[[x]]]/2;
        tauout <- rgamma(1, shape = a + length(storseq[[x]])/2, rate = apart + b);
        return(tauout)})
      tau <- c(NA,tau);
      matlist <- unlist(c(B_inv, c(sapply(2:length(dim_q),function(x) rep(tau[x],dim_q[x])))))
      Big_Prec <- diag(matlist, nrow = length(matlist))
    }
    
    # Save posterior sample 
    if(samp>warmup){
      for(k in 1:l_q){
        beta.samps[[k]][samp-warmup,] <- beta_vec[storseq[[k]]]
        if(k > 1){
          tau.samps[[k-1]][samp-warmup] <- tau[k]
        }}
      delta.samps <- deltas + delta.samps
      epsilon.samps[samp-warmup] <- epsilons
      gamma.samps[samp-warmup] <- gammas
    } 
  }
  return(list(Beta_Samps = beta.samps,Tau_samps=tau.samps, Delta_Samps = delta.samps/(nsamps-warmup),
              Epsilon_Samps = epsilon.samps)) 
}

P_Eps_Eval_Simp <- function(a_T, b_T, length_Q,deltas){
  rbeta(1,a_T + sum(deltas), b_T + (length_Q - sum(deltas)))
}

# provides estimates for parameters yielding a suitable 
# beta distribution to reflect the rate of unaccounted for 
# vaccinations 
Estimate_Distribution_Unaccount_Vaccines <- function(mode, low, up, 
                                                     q_low, q_up){
  try_qtbeta <- function(alpha, mode, low, up,q_low, q_up){
    beta <- (alpha*(1-mode) + (2*mode -1))/mode
    qs <- c(q_low[1], q_up[1])
    percents <- c(q_low[2], q_up[2])
    ptbeta(qs[2], alpha = alpha, beta = beta, a = low, b = up)-ptbeta(qs[1], alpha = alpha, beta = beta, a = low, b = up)-diff(percents)
  }
  
  alphaest <- uniroot(try_qtbeta, interval= c(1,1000),
                      q_low = q_low, q_up = q_up, mode= mode, low = low, up = up, extendInt = "upX")$root
  
  beta_est <- (alphaest*(1-mode) + (2*mode -1))/mode
  
  return(c(alphaest, beta_est))
}


# Data_Frames: list of X - fixed effects matrix, A_1 - random effect matrix, ..., A_R - Rth random effect matrix 
# Y = number of sites length vector of vaccination counts per site
# N = number of sites length vector with total number of vaccine-eligible veterans at each site
# nsamps = number of posterior samples 
# warmup = burnin stage (discarded posterior samples, enough for the sampler to 'reach convergence')
# hyp = vector of prior rate and scale parameters for inverse gamma prior on random effect variances
# these two numbers are used for all hiearchical levels/random effects' variances
# Beta_Params: a Vector of Length 7 with the following properties: 
# (Most likely value for unaccounted for vaccine rate, lower bound on unaccounted for vaccine rate, 
# upper bound on unaccounted for vaccine rate, lower %-ile value, lower %-ile %, upper %-ile value, upper %-ile %)
# this prior is to overcome how suggestive the prior distribution was in the method 'site_gamma1'; 
# the normal distribution is much too suggestive and not very clear on how to specify the variance
# going through the rate of unaccounted for vaccinations seems to be more applicable and could be drawn 
# easier from admin/US level data. 
# also, Beta_Parameter can be a matrix of (# of sites as row count by 7 columns) all hypothesized 
# unaccounted for vaccine rate at each site. 
# underreport = TRUE indicates you are using the proposed method, underreport = FALSE indicates you are simply
# using hierarchical logistic regression in a bayesian framework, no imputation step 
# VA_Only = FALSE indicates the focus is on vaccination rates in general
# VA_Only = TRUE indicates the focus is on vaccination rates in the VA only. 

site_pgamma2 = function(Data_Frames,Y,N,nsamps=4000,warmup=500, hyp = c(0.01, 0.01),underreport = TRUE, VA_only = FALSE,
                        Beta_Params = c(0.05, 0.02, 0.20, 0.03, 0.05, 0.15, 0.95)){
  # this is on the site level, so nrow(X) = number of sites 
  # fixed effects first in Data_Frames list
  a = hyp[1]; b = hyp[2]
  
  # estimating parameters for beta distribution for draw
  # of vaccination rate thats unaccounted for
  if(!is.matrix(Beta_Params) | is.null(nrow(Beta_Params))){
    mode = Beta_Params[1]; low = Beta_Params[2]; up = Beta_Params[3]; 
    q_low = Beta_Params[4:5]; q_up = Beta_Params[6:7]
Beta_Samp_Params <- Estimate_Distribution_Unaccount_Vaccines(mode, low, up, 
                                                         q_low, q_up)
Beta_Samp_Params <- as.matrix(t(Beta_Samp_Params), nrow =1)
  }else{
    mode = Beta_Params[,1]; low = Beta_Params[,2]; up = Beta_Params[,3]; 
    q_low = Beta_Params[,4:5]; q_up = Beta_Params[,6:7]
Beta_Samp_Params <- t(sapply(1:nrow(Beta_Params), function(x) Estimate_Distribution_Unaccount_Vaccines(mode[x], low[x], up[x], 
                                                              q_low[x,], q_up[x,])))
  }


  l_q <- length(Data_Frames)
  dim_q <- sapply(1:l_q, function(x) ncol(Data_Frames[[x]]))
  
  j = 1; storseq <- list()
  
  for(k in 1:l_q){
    dim1 <- dim_q[k]
    if(k == 1){seq1 <- j:dim1}else{
      seq1 <- j:sum(dim_q[1:k])
    }
    storseq[[k]] <- seq1
    j <- max(seq1)+1
  }
  
  tot_effs <- sum(dim_q)
  n_site <- W <- nrow(Data_Frames[[1]])
  X <- do.call(cbind,Data_Frames)
  
  # random jitter in outcomes:
  if(underreport){
    if(VA_only){
      delta <- ceiling((N-Y-1)*rtbeta(W, alpha = Beta_Samp_Params[,1], beta = Beta_Samp_Params[,2], a = low, b = up))
      N_Star <- N-delta # what we think if VA site only
      Y_Star = Y
    }else{
      delta <- ceiling((N-Y-1)*rtbeta(W, alpha = Beta_Samp_Params[,1], beta = Beta_Samp_Params[,2], a = low, b = up))
      N_Star <- N
      Y_Star = Y + delta # what we think if we include all possible data sources
    }
  }else{
    Y_Star = Y; N_Star = N
  }
  kappa = (Y_Star-0.5*N_Star)
  # Storage object
  beta.samps <- lapply(1:l_q, function(y) matrix(NA,nrow=nsamps,ncol = dim_q[y]))
  if(l_q > 1){
    tau.samps <- lapply(1:(l_q-1), function(x) rep(NA,nsamps))
    tau = rep(0.5, l_q) 
    tau = tau[-1]
  }else{
    tau.samps <- NULL; 
    tau = NULL
  }
  
  #inits
  beta_vec <- rep(1, tot_effs)
  B_inv <- rep(0.01, dim_q[1])
  matlist <- unlist(c(B_inv, c(sapply(2:length(dim_q),function(x) rep(tau[x-1],dim_q[x])))))
  Big_Prec <- diag(matlist, nrow = length(matlist))
  
  # Gibbs Sampler  
  for(samp in 1:(nsamps+warmup)){
    # Sample alpha
    #print(paste0("Posterior run ", samp, ' started'))
    
    omega = rpg.devroye(n_site,N_Star,X%*%beta_vec)        # Update auxillary parameters
    
    V_Beta = chol2inv(chol(crossprod(X,(X*omega)) + Big_Prec))
    m_beta = as.matrix(V_Beta%*%crossprod(X,kappa))
    # Mean vector
    beta_vec = (m_beta + t(chol(V_Beta))%*%rnorm(ncol(X)))# Unconstrained
    
    if(l_q > 1){
      #sample from posterior for tau
      tau = sapply(2:l_q, function(x){
        apart <- t(beta_vec[storseq[[x]]])%*%beta_vec[storseq[[x]]]/2;
        tauout <- rgamma(1, shape = a + length(storseq[[x]])/2, rate = apart + b);
        return(tauout)})
      tau <- c(NA,tau);
      matlist <- unlist(c(B_inv, c(sapply(2:length(dim_q),function(x) rep(tau[x],dim_q[x])))))
      Big_Prec <- diag(matlist, nrow = length(matlist))
    }
    
    # random jitter in outcomes:
    if(underreport){
      if(VA_only){
        delta <- ceiling((N-Y-1)*rtbeta(W, alpha = Beta_Samp_Params[,1], beta = Beta_Samp_Params[,2], a = low, b = up))
        N_Star <- N-delta # what we think if VA site only
        Y_Star = Y
      }else{
        delta <- ceiling((N-Y-1)*rtbeta(W, alpha = Beta_Samp_Params[,1], beta = Beta_Samp_Params[,2], a = low, b = up))
        N_Star <- N
        Y_Star = Y + delta # what we think if we include all possible data sources
      }
      # we need to restrict this fraction to be smaller than 1. 
      # this is how. 
      kappa = (Y_Star-0.5*N_Star)
    }
    # Save posterior sample 
    if(samp>warmup){
      for(k in 1:l_q){
        beta.samps[[k]][samp-warmup,] <- beta_vec[storseq[[k]]]
        if(k > 1){
          tau.samps[[k-1]][samp-warmup] <- tau[k]
        }}
    } 
  }
  return(list(Beta_Samps = beta.samps,Tau_samps=tau.samps)) 
}

# Data_Frames: list of X - fixed effects matrix, A_1 - random effect matrix, ..., A_R - Rth random effect matrix 
# Y = number of sites length vector of vaccination counts per site
# N = number of sites length vector with total number of vaccine-eligible veterans at each site
# nsamps = number of posterior samples 
# warmup = burnin stage (discarded posterior samples, enough for the sampler to 'reach convergence')
# hyp = vector of prior rate and scale parameters for inverse gamma prior on random effect variances
# these two numbers are used for all hiearchical levels/random effects' variances
# Beta_Params: a Vector of Length 7 with the following properties: 
# (Most likely value for unaccounted for vaccine rate, lower bound on unaccounted for vaccine rate, 
# upper bound on unaccounted for vaccine rate, lower %-ile value, lower %-ile %, upper %-ile value, upper %-ile %)
# this prior is to overcome how suggestive the prior distribution was in the method 'site_gamma1'; 
# the normal distribution is much too suggestive and not very clear on how to specify the variance
# going through the rate of unaccounted for vaccinations seems to be more applicable and could be drawn 
# easier from admin/US level data. 
# also, Beta_Parameter can be a matrix of (# of sites as row count by 7 columns) all hypothesized 
# unaccounted for vaccine rate at each site. 
# underreport = TRUE indicates you are using the proposed method, underreport = FALSE indicates you are simply
# using hierarchical logistic regression in a bayesian framework, no imputation step 
# VA_Only = FALSE indicates the focus is on vaccination rates in general
# VA_Only = TRUE indicates the focus is on vaccination rates in the VA only. 
# if Beta_Params = NULL then no

site_pgamma3 = function(Data_Frames,Y,N,nsamps=4000,warmup=500, hyp = c(0.01, 0.01),underreport = TRUE,
                        Beta_Params = c(0.05, 0.02, 0.20, 0.03, 0.05, 0.15, 0.95),
                        Csi_Params = c(0.07, 0.01, 0.25, 0.03, 0.05, 0.20, 0.95)){
  # this is on the site level, so nrow(X) = number of sites 
  # fixed effects first in Data_Frames list
  a = hyp[1]; b = hyp[2]
  flag_beta = F; flag_beta2 = F
  # estimating parameters for beta distribution for draw
  # of vaccination rate thats unaccounted for
  
  if(underreport == TRUE){
  if(!is.null(Beta_Params)){
  if(!is.matrix(Beta_Params) | is.null(nrow(Beta_Params))){
    mode = Beta_Params[1]; low = Beta_Params[2]; up = Beta_Params[3]; 
    q_low = Beta_Params[4:5]; q_up = Beta_Params[6:7]
    Beta_Samp_Params <- Estimate_Distribution_Unaccount_Vaccines(mode, low, up, 
                                                                 q_low, q_up)
    Beta_Samp_Params <- as.matrix(t(Beta_Samp_Params), nrow =1)
  }else{
    mode = Beta_Params[,1]; low = Beta_Params[,2]; up = Beta_Params[,3]; 
    q_low = Beta_Params[,4:5]; q_up = Beta_Params[,6:7]
    Beta_Samp_Params <- t(sapply(1:nrow(Beta_Params), function(x) Estimate_Distribution_Unaccount_Vaccines(mode[x], low[x], up[x],q_low[x,], q_up[x,])))
  }
  }else{flag_beta = T}
  
  if(!is.null(Csi_Params)){
    if(!is.matrix(Csi_Params) | is.null(nrow(Csi_Params))){
      mode = Csi_Params[1]; low2 = Csi_Params[2]; up2 = Csi_Params[3]; 
      q_low = Csi_Params[4:5]; q_up = Csi_Params[6:7]
      Beta_Samp_Params2 <- Estimate_Distribution_Unaccount_Vaccines(mode, low2, up2, 
                                                                   q_low, q_up)
      Beta_Samp_Params2 <- as.matrix(t(Beta_Samp_Params2), nrow =1)
    }else{
      mode = Csi_Params[,1]; low2 = Csi_Params[,2]; up2 = Csi_Params[,3]; 
      q_low = Csi_Params[,4:5]; q_up = Csi_Params[,6:7]
      Beta_Samp_Params2 <- t(sapply(1:nrow(Csi_Params), function(x) Estimate_Distribution_Unaccount_Vaccines(mode[x], low2[x], up2[x],q_low[x,], q_up[x,])))
    }
  }else{flag_beta2 = T}
}
  
  l_q <- length(Data_Frames)
  dim_q <- sapply(1:l_q, function(x) ncol(Data_Frames[[x]]))
  
  j = 1; storseq <- list()
  
  for(k in 1:l_q){
    dim1 <- dim_q[k]
    if(k == 1){seq1 <- j:dim1}else{
      seq1 <- j:sum(dim_q[1:k])
    }
    storseq[[k]] <- seq1
    j <- max(seq1)+1
  }
  
  tot_effs <- sum(dim_q)
  n_site <- W <- nrow(Data_Frames[[1]])
  X <- do.call(cbind,Data_Frames)
  
  # random jitter in outcomes:
  if(underreport){
      if(!flag_beta2){
      delta_2 <- ceiling((N-Y-1)*rtbeta(W, alpha = Beta_Samp_Params2[,1],beta = Beta_Samp_Params2[,2], a = low2, b = up2))
      }else{delta_2 = 0}
      N_Star <- N-delta_2 # at baseline, who are misclassified as not vaccinated but vaccinated at or before baseline
      if(!flag_beta){
        delta <- ceiling((N_Star-Y-1)*rtbeta(W, alpha = Beta_Samp_Params[,1], beta = Beta_Samp_Params[,2], a = low, b = up))
      }else{delta = 0}
      # changed to N-delta-delta_2 on august 31st 2023 because i forgot we are also moving delta_2 to numerator FRom denominator!
      Y_Star = Y + delta # during study period, those not flagged as vaccinated 
  }else{
    Y_Star = Y; N_Star = N
  }
  
  kappa = (Y_Star-0.5*N_Star)
  # Storage object
  beta.samps <- lapply(1:l_q, function(y) matrix(NA,nrow=nsamps,ncol = dim_q[y]))
  if(l_q > 1){
    tau.samps <- lapply(1:(l_q-1), function(x) rep(NA,nsamps))
    tau = rep(0.5, l_q) 
    tau = tau[-1]
  }else{
    tau.samps <- NULL; 
    tau = NULL
  }
  
  #inits
  beta_vec <- rep(1, tot_effs)
  B_inv <- rep(1, dim_q[1])
  matlist <- unlist(c(B_inv, c(sapply(2:length(dim_q),function(x) rep(tau[x-1],dim_q[x])))))
  Big_Prec <- diag(matlist, nrow = length(matlist))
  
  # Gibbs Sampler  
  for(samp in 1:(nsamps+warmup)){
    # Sample alpha
    #print(paste0("Posterior run ", samp, ' started'))
    
    omega = rpg.devroye(n_site,N_Star,X%*%beta_vec)        # Update auxillary parameters
    
    V_Beta = chol2inv(chol(crossprod(X,(X*omega)) + Big_Prec))
    m_beta = as.matrix(V_Beta%*%crossprod(X,kappa))
    # Mean vector
    beta_vec = (m_beta + t(chol(V_Beta))%*%rnorm(ncol(X)))# Unconstrained
    
    if(l_q > 1){
      #sample from posterior for tau
      tau = sapply(2:l_q, function(x){
        apart <- t(beta_vec[storseq[[x]]])%*%beta_vec[storseq[[x]]]/2;
        tauout <- rgamma(1, shape = a + length(storseq[[x]])/2, rate = apart + b);
        return(tauout)})
      tau <- c(NA,tau);
      matlist <- unlist(c(B_inv, c(sapply(2:length(dim_q),function(x) rep(tau[x],dim_q[x])))))
      Big_Prec <- diag(matlist, nrow = length(matlist))
    }
    
    # random jitter in outcomes:
    if(underreport){
      if(!flag_beta2){
        delta_2 <- ceiling((N-Y-1)*rtbeta(W, alpha = Beta_Samp_Params2[,1],beta = Beta_Samp_Params2[,2], a = low2, b = up2))
      }else{delta_2 = 0}
      N_Star <- N-delta_2 # at baseline, who are misclassified as not vaccinated but vaccinated at or before baseline
      if(!flag_beta){
        delta <- ceiling((N_Star-Y-1)*rtbeta(W, alpha = Beta_Samp_Params[,1], beta = Beta_Samp_Params[,2], a = low, b = up))
      }else{delta = 0}
      # changed to be sequental deduction in N_Star, then take additional Y from N_Star
      # because that's what we are assuming. 
      Y_Star = Y + delta # during study period, those not flagged as vaccinated 
    }else{
      Y_Star = Y; N_Star = N
    }
    
    kappa = (Y_Star-0.5*N_Star)
    # Save posterior sample 
    if(samp>warmup){
      for(k in 1:l_q){
        beta.samps[[k]][samp-warmup,] <- beta_vec[storseq[[k]]]
        if(k > 1){
          tau.samps[[k-1]][samp-warmup] <- tau[k]
        }}
    } 
  }
  return(list(Beta_Samps = beta.samps,Tau_samps=tau.samps)) 
}

rho1_plot <- function(beta0, beta1, rho1s, rho2s, rho2){
p0 <- (1+exp(-(beta0+ beta1)))^(-1)
numerator = p0*((1+rho2s)*(1+rho2s*(1-(1-rho1s)*p0)) - rho2*(1-rho1s)*p0) - p0*(1-rho1s)
denom = (1+rho2)*(1 + rho2s*(1 - (1-rho1s)*p0)) - (1+rho2)*(p0*(1-rho1s))
outd <- numerator/denom
return(outd)
}


logit <- function(x){log(x)-log(1-x)}
expit <- function(x){exp(x)/(1+exp(x))}

Grab_ConvertP0True <- function(p0o, rho1s, rho2s){
  try_p0 <- function(p0,p0o, rho1s, rho2s){
  (((1-rho1s)*p0)/(1+rho2s*(1-p0)))- p0o
  }
  
  alphaest <- uniroot(try_p0, interval=c(0, 1),
                      p0o = p0o, rho1s = rho1s, rho2s = rho2s, extendInt = 'no')$root
  return(alphaest)
}

Grab_ConvertP0True2 <- function(p0o, rho1s, rho2s){
  try_p0 <- function(p0,p0o, rho1s, rho2s){
    x1 = 1 - (p0o/(1-rho2s*(1-p0o)))
    x2 = (1-p0o)/(1-rho2s*(1-p0o))
    ((p0 - rho1s*(x1))/(1 + rho2s*(x2)))- p0o
  }
  
  alphaest <- uniroot(try_p0, interval=c(0, 1),
                      p0o = p0o, rho1s = rho1s, rho2s = rho2s, extendInt = 'no')$root
  return(alphaest)
}

Recover_Observed_P0 <- function(p0, rho1, rho2){
  try_p0 <- function(p0o,p0, rho1s, rho2s){
    x1 = 1 - (p0o/(1-rho2*(1-p0o)))
    x2 = (1-p0o)/(1-rho2*(1-p0o))
    ((p0 - rho1*(x1))/(1 + rho2*(x2)))- p0o
  }
  
  alphaest <- uniroot(try_p0, interval=c(0, 1),
                      p0 = p0, rho1 = rho1, rho2 = rho2, extendInt = 'no')$root
  return(alphaest)
}
# functions from TruncatedDistribution Package 

rtbeta <- function(n, alpha, beta, a=0, b=1)
{
  stopifnot(n > 0 & all(beta > 0) & all(alpha > 0))
  x <- runif(n)
  Fa <- pbeta(a, alpha, beta)
  Fb <- pbeta(b, alpha, beta)
  y <- (1-x)*Fa + x*Fb
  return(qbeta(y, alpha, beta))
}

dtbeta <- function(x, alpha, beta, a=0, b=1)
{
  stopifnot( all(alpha > 0) & all(beta > 0) )
  Fa <- pbeta(a, alpha, beta)
  Fb <- pbeta(b, alpha, beta)
  y <- dbeta(x, alpha, beta)
  inda <- which(x < a)
  indb <- which(x > b)
  if(length(inda) > 0) y[inda] <- 0
  if(length(indb) > 0) y[indb] <- 0
  return(y/(Fb-Fa))
}

ptbeta <- function(q, alpha, beta, a=0, b=1)
{
  stopifnot( all( alpha > 0 ) & all(beta > 0) )
  Fa <- pbeta(a, alpha, beta)
  Fb <- pbeta(b, alpha, beta)
  p <- ( pbeta(q, alpha, beta) - Fa ) / ( Fb - Fa )
  inda <- which(q < a)
  indb <- which(q > b)
  if(length(inda) > 0) p[inda] <- 0
  if(length(indb) > 0) p[indb] <- 1
  return(p)
}

qtbeta <- function(p, alpha, beta, a=0, b=1)
{
  stopifnot( all(p >= 0 & p <= 1) &  all( alpha > 0 ) & all(beta > 0) )
  Fa <- pbeta(a, alpha, beta)
  Fb <- pbeta(b, alpha, beta)
  pNew <- p * (Fb - Fa) + Fa
  x <- qbeta(pNew, alpha, beta)
  return(x)
}

# adding in the change for eligibility (4/22/2024)
# Data_Frames: list of X - fixed effects matrix, A_1 - random effect matrix, ..., A_R - Rth random effect matrix 
# Y = number of sites length vector of vaccination counts per site
# N = number of sites length vector with total number of vaccine-eligible veterans at each site
# nsamps = number of posterior samples 
# warmup = burnin stage (discarded posterior samples, enough for the sampler to 'reach convergence')
# hyp = vector of prior rate and scale parameters for inverse gamma prior on random effect variances
# these two numbers are used for all hiearchical levels/random effects' variances
# Beta_Params: a Vector of Length 7 with the following properties: 
# (Most likely value for unaccounted for vaccine rate, lower bound on unaccounted for vaccine rate, 
# upper bound on unaccounted for vaccine rate, lower %-ile value, lower %-ile %, upper %-ile value, upper %-ile %)
# this prior is to overcome how suggestive the prior distribution was in the method 'site_gamma1'; 
# the normal distribution is much too suggestive and not very clear on how to specify the variance
# going through the rate of unaccounted for vaccinations seems to be more applicable and could be drawn 
# easier from admin/US level data. 
# also, Beta_Parameter can be a matrix of (# of sites as row count by 7 columns) all hypothesized 
# unaccounted for vaccine rate at each site. 
# underreport = TRUE indicates you are using the proposed method, underreport = FALSE indicates you are simply
# using hierarchical logistic regression in a bayesian framework, no imputation step 
# VA_Only = FALSE indicates the focus is on vaccination rates in general
# VA_Only = TRUE indicates the focus is on vaccination rates in the VA only. 
# if Beta_Params = NULL then no

site_pgamma4 = function(Data_Frames,Y,N,nsamps=4000,warmup=500, hyp = c(0.01, 0.01),underreport = TRUE,
                        Beta_Params = c(0.05, 0.02, 0.20, 0.03, 0.05, 0.15, 0.95),
                        Csi_Params = c(0.07, 0.01, 0.25, 0.03, 0.05, 0.20, 0.95),
                        Csi_Params2 = c(0.07, 0.01, 0.25, 0.03, 0.05, 0.20, 0.95)){
  # this is on the site level, so nrow(X) = number of sites 
  # fixed effects first in Data_Frames list
  a = hyp[1]; b = hyp[2]
  flag_beta = F; flag_beta2 = F; flag_beta3 = F
  # estimating parameters for beta distribution for draw
  # of vaccination rate thats unaccounted for
  
  if(underreport == TRUE){
    if(!is.null(Beta_Params)){
      if(!is.matrix(Beta_Params) | is.null(nrow(Beta_Params))){
        mode = Beta_Params[1]; low = Beta_Params[2]; up = Beta_Params[3]; 
        q_low = Beta_Params[4:5]; q_up = Beta_Params[6:7]
        Beta_Samp_Params <- Estimate_Distribution_Unaccount_Vaccines(mode, low, up, 
                                                                     q_low, q_up)
        Beta_Samp_Params <- as.matrix(t(Beta_Samp_Params), nrow =1)
      }else{
        mode = Beta_Params[,1]; low = Beta_Params[,2]; up = Beta_Params[,3]; 
        q_low = Beta_Params[,4:5]; q_up = Beta_Params[,6:7]
        Beta_Samp_Params <- t(sapply(1:nrow(Beta_Params), function(x) Estimate_Distribution_Unaccount_Vaccines(mode[x], low[x], up[x],q_low[x,], q_up[x,])))
      }
    }else{flag_beta = T}
    
    if(!is.null(Csi_Params)){
      if(!is.matrix(Csi_Params) | is.null(nrow(Csi_Params))){
        mode = Csi_Params[1]; low2 = Csi_Params[2]; up2 = Csi_Params[3]; 
        q_low = Csi_Params[4:5]; q_up = Csi_Params[6:7]
        Beta_Samp_Params2 <- Estimate_Distribution_Unaccount_Vaccines(mode, low2, up2, 
                                                                      q_low, q_up)
        Beta_Samp_Params2 <- as.matrix(t(Beta_Samp_Params2), nrow =1)
      }else{
        mode = Csi_Params[,1]; low2 = Csi_Params[,2]; up2 = Csi_Params[,3]; 
        q_low = Csi_Params[,4:5]; q_up = Csi_Params[,6:7]
        Beta_Samp_Params2 <- t(sapply(1:nrow(Csi_Params), function(x) Estimate_Distribution_Unaccount_Vaccines(mode[x], low2[x], up2[x],q_low[x,], q_up[x,])))
      }
    }else{flag_beta2 = T}
    
    if(!is.null(Csi_Params2)){
      if(!is.matrix(Csi_Params2) | is.null(nrow(Csi_Params2))){
        mode = Csi_Params2[1]; low3 = Csi_Params2[2]; up3 = Csi_Params2[3]; 
        q_low = Csi_Params2[4:5]; q_up = Csi_Params2[6:7]
        Beta_Samp_Params3 <- Estimate_Distribution_Unaccount_Vaccines(mode, low3, up3, 
                                                                      q_low, q_up)
        Beta_Samp_Params3 <- as.matrix(t(Beta_Samp_Params3), nrow =1)
      }else{
        mode = Csi_Params2[,1]; low3 = Csi_Params2[,2]; up3 = Csi_Params2[,3]; 
        q_low = Csi_Params2[,4:5]; q_up = Csi_Params2[,6:7]
        Beta_Samp_Params3 <- t(sapply(1:nrow(Csi_Params2), function(x) Estimate_Distribution_Unaccount_Vaccines(mode[x], low3[x], up3[x],q_low[x,], q_up[x,])))
      }
    }else{flag_beta3 = T}
  }
  
  l_q <- length(Data_Frames)
  dim_q <- sapply(1:l_q, function(x) ncol(Data_Frames[[x]]))
  
  j = 1; storseq <- list()
  
  for(k in 1:l_q){
    dim1 <- dim_q[k]
    if(k == 1){seq1 <- j:dim1}else{
      seq1 <- j:sum(dim_q[1:k])
    }
    storseq[[k]] <- seq1
    j <- max(seq1)+1
  }
  
  tot_effs <- sum(dim_q)
  n_site <- W <- nrow(Data_Frames[[1]])
  X <- do.call(cbind,Data_Frames)
  
  # random jitter in outcomes:
  if(underreport){
    if(!flag_beta3){
      delta_3 <- ceiling((Y)*rtbeta(W, alpha = Beta_Samp_Params3[,1],beta = Beta_Samp_Params3[,2], a = low3, b = up3))
      Z = Y-delta_3
      }else{delta_3 = 0; Z = Y} 

    if(!flag_beta2){
      delta_2 <- ceiling((N-Z)*rtbeta(W, alpha = Beta_Samp_Params2[,1],beta = Beta_Samp_Params2[,2], a = low2, b = up2))
    }else{delta_2 = 0}
    N_Star <- N-delta_2-delta_3 # at baseline, who are misclassified as not vaccinated but vaccinated at or before baseline
    # this is all eligibility corrections above
    if(!flag_beta){
      delta <- ceiling((N_Star-Z)*rtbeta(W, alpha = Beta_Samp_Params[,1], beta = Beta_Samp_Params[,2], a = low, b = up))
    }else{delta = 0}
    # changed to N-delta-delta_2 on august 31st 2023 because i forgot we are also moving delta_2 to numerator FRom denominator!
    Y_Star = Z + delta # during study period, those not flagged as vaccinated 
  }else{
    Y_Star = Y; N_Star = N
  }

  
  kappa = (Y_Star-0.5*N_Star)
  # Storage object
  beta.samps <- lapply(1:l_q, function(y) matrix(NA,nrow=nsamps,ncol = dim_q[y]))
  if(l_q > 1){
    tau.samps <- lapply(1:(l_q-1), function(x) rep(NA,nsamps))
    tau = rep(0.5, l_q) 
    tau = tau[-1]
  }else{
    tau.samps <- NULL; 
    tau = NULL
  }
  
  #inits
  beta_vec <- rep(1, tot_effs)
  B_inv <- rep(1, dim_q[1])
  matlist <- unlist(c(B_inv, c(sapply(2:length(dim_q),function(x) rep(tau[x-1],dim_q[x])))))
  Big_Prec <- diag(matlist, nrow = length(matlist))
  
  # Gibbs Sampler  
  for(samp in 1:(nsamps+warmup)){
    # Sample alpha
    #print(paste0("Posterior run ", samp, ' started'))
    
    omega = rpg.devroye(n_site,N_Star,X%*%beta_vec)# Update auxillary parameters
    
    V_Beta = chol2inv(chol(crossprod(X,(X*omega)) + Big_Prec))
    m_beta = as.matrix(V_Beta%*%crossprod(X,kappa))
    # Mean vector
    beta_vec = (m_beta + t(chol(V_Beta))%*%rnorm(ncol(X)))# Unconstrained
    
    if(l_q > 1){
      #sample from posterior for tau
      tau = sapply(2:l_q, function(x){
        apart <- t(beta_vec[storseq[[x]]])%*%beta_vec[storseq[[x]]]/2;
        tauout <- rgamma(1, shape = a + length(storseq[[x]])/2, rate = apart + b);
        return(tauout)})
      tau <- c(NA,tau);
      matlist <- unlist(c(B_inv, c(sapply(2:length(dim_q),function(x) rep(tau[x],dim_q[x])))))
      Big_Prec <- diag(matlist, nrow = length(matlist))
    }
    
    # random jitter in outcomes:
    if(underreport){
      if(!flag_beta3){
        delta_3 <- ceiling((Y)*rtbeta(W, alpha = Beta_Samp_Params3[,1],beta = Beta_Samp_Params3[,2], a = low3, b = up3))
        Z = Y-delta_3
      }else{delta_3 = 0; Z = Y} 
    if(!flag_beta2){
      delta_2 <- ceiling((N-Z)*rtbeta(W, alpha = Beta_Samp_Params2[,1],beta = Beta_Samp_Params2[,2], a = low2, b = up2))
    }else{delta_2 = 0}
    N_Star <- N-delta_2-delta_3 # at baseline, who are misclassified as not vaccinated but vaccinated at or before baseline
    # this is all eligibility corrections above
    if(!flag_beta){
      delta <- ceiling((N_Star-Z)*rtbeta(W, alpha = Beta_Samp_Params[,1], beta = Beta_Samp_Params[,2], a = low, b = up))
    }else{delta = 0}
    # changed to N-delta-delta_2 on august 31st 2023 because i forgot we are also moving delta_2 to numerator FRom denominator!
    Y_Star = Z + delta # during study period, those not flagged as vaccinated 
  }else{
    Y_Star = Y; N_Star = N
  }
    kappa = (Y_Star-0.5*N_Star)
    
    # Save posterior sample 
    if(samp>warmup){
      for(k in 1:l_q){
        beta.samps[[k]][samp-warmup,] <- beta_vec[storseq[[k]]]
        if(k > 1){
          tau.samps[[k-1]][samp-warmup] <- tau[k]
        }}
    } 
  }
  return(list(Beta_Samps = beta.samps,Tau_samps=tau.samps)) 
}
