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

logit <- function(x){log(x)-log(1-x)}
expit <- function(x){exp(x)/(1+exp(x))}

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
# Csi_Params: a vector of Length 7 with the same properties as Beta_Params but for the misclassified successes 
# in terms of eligibility 
# Csi_Params2: same as Csi_Params however for the misclassified failures in terms of eligibility
# underreport = TRUE indicates you are using the proposed method, underreport = FALSE indicates you are simply
# using hierarchical logistic regression in a bayesian framework, no imputation step 
# if Beta_Params = NULL then no misclassification in numerator
# same for Csi_Params and Csi_Params2

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
