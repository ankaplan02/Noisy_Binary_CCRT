# more sims

# simulations based on given observed rate and rates of misclassification 
# then give true value, i.e., implying the true rate was this. 
rawdat <- readRDS("rawdat.rds")
source("FuzzyBinomialDataPGammaModel.R")

##################### constructing the appropriate data frames from rawdat.####################
################# matrices X, A1, A2, A_....#####################################

ngroup1 <- length(table(rawdat$Group1)) # 10 health care systems
ngroup2 <- length(table(rawdat$Group2)) # 90 clinics
n_group_per1 <- table(rawdat$Group1) # count of clinics per HCS

X <- cbind(rep(1, nrow(rawdat)), as.numeric(rawdat$StudyArm))

N_obs <- ceiling(rawdat$Frequency)

n_site <- nrow(rawdat)

id1 <- unlist(lapply(1:length(n_group_per1), function(x) rep(x, n_group_per1[x])))
id2 <- 1:nrow(rawdat)

A <- diag(1, nrow = nrow(rawdat)); 
A2 <- model.matrix(~-1+factor(rawdat$Group1))

Data_Frames = list(as.matrix(X), A,A2)


### misclassification rates specified 

rho2 = .04; rho1 = .07; # 4% in general have misclassified eligibility;
# 7% in control arm failures should be vaccinated and eligible
rho22 = .04; rho12 = .13 
# rho1 = 13% of eligible failures should be vaccinated (Y=1s) in treatment arm

# in truth

# this will be the sample size per site (N_star)
# a vector; we'll have to do this with a rho1/rho2 column vector
# rho1/rho2 at places where an arm is a treatment arm or control arm
require(dplyr)
# possibly do this in a DPLYR step where ifelse statements if treatment arm or not
ratep = c(.30, .33) # observed vaccination rate 
# observed vaccination rate AMONG control group and vaccination arm 

# observed odds ratio is 
obsOR <- ratep[2]*(1-ratep[1])/(ratep[1]*(1-ratep[2]))
# true odds ratio 
# based on eligible percentages (so given we are in the eligible pool)

# on average the rates in truth
new_ratec <- ratep[1] + (1-ratep[1])*rho1
new_ratet <- ratep[2] + (1-ratep[2])*rho12
ratep2 <- c(new_ratec, new_ratet)
# true odds ratio is 
#
trueOR <- ratep2[2]*(1-ratep2[1])/(ratep2[1]*(1-ratep2[2]))

# specifying the misclassification rates for each arm. We can do this via a list of matrices. We will 
# specify a rho specification for each site

# 13% among eligible failures in treatment arm should be relabeled as successes
# 7% of eligible failures in control arm should be relabeled as successes
# both arms have 4 percent eligibility N are incorrect

rho1set <- c(.13, 0.12, 0.14) # specified to be treatment rhos in truth (3 options)
bnds1 <- cbind(rho1set, c(.11, .11, .11), c(0.15,.15,.15), 
               rho1set-.005,rep(0.05, 3), 
               rho1set + .005, rep(0.95, 3)) # this cbind creates the specifications rho-hat, lowest and highest
# assumed values, with a 95% certainty around rho-hat by 1% density (0.005 + 0.005)
bnds12 <- c(0.07, .06, 0.08, 0.065, 0.05, 0.075, 0.95) # for control its specifications
# (rho-hat, lowest, highest, 5%-ile of rho, 5%, 95%-ile of rho, 95%)
# changing the percentiles to not be quite symmetric around rho-hat yields left-skewed/right-skewed distributions
# REALLY TIGHT AROUND THEIR true values
csi1_and2 <- c(0.04, .03, 0.05, 0.035, 0.05, 0.045, 0.95) # this will be used for both arms; the eligibility misclassification
# eligibility will effect both Y and N-Y people in each arm equally 

##################### Specifying differential misclassification between arms, for each site#################### 

Beta_Set <- list()
Beta_Set[[1]] <- matrix(0, nrow = 7, ncol = nrow(X))
Beta_Set[[2]] <- matrix(0, nrow = 7, ncol = nrow(X))
Beta_Set[[3]] <- matrix(0, nrow = 7, ncol = nrow(X))

# the control arm is the one changing in each set of sims
Beta_Set[[1]][,which(X[,2]==0)] <- bnds12 # we assume all control clinics are modeled by the same misclassification rate in each scenario
Beta_Set[[1]][,which(X[,2]==1)] <- bnds1[1,] # when X[,2] == 1, these clinics are treatment intervention assigned clinics
Beta_Set[[2]][,which(X[,2]==0)] <- bnds12
Beta_Set[[2]][,which(X[,2]==1)] <- bnds1[2,]
Beta_Set[[3]][,which(X[,2]==0)] <- bnds12
Beta_Set[[3]][,which(X[,2]==1)] <- bnds1[3,]

Beta_Set <- lapply(Beta_Set, function(x) t(x)) # transpose each matrix. 

# each element of Beta_Set list is a misclassification we assume for modeling, not what we assume in data generation

# three sets

n_set_beta <- length(Beta_Set)
# current three sims have this set up: 

ICC_site = 0.07; # intra-class correlation values used
ICC_hcl <- 0.02;
# works for one to two ICCs (two to three levels)

sigs <- convertICCtoVar(c(ICC_site, ICC_hcl)) # convert ICC to variance values (works only for 3-level intercept or 2-level intercept)
sigma_site <- sigs[1]; sigma_hcl <- sigs[2] # random effect variances for multi-level hierarchy

#############################################################################################################
########################### Generating Outcome Data############################################################

Data <- rawdat %>% mutate(TRTArm = X[,2]) %>% rowwise() %>%
  mutate(N_true = rbinom(1, ceiling(Frequency), 1-rho2), # generate true N (eligible patients)
         Rhos = ifelse(TRTArm == 1, rho12, rho1), # assign the correct Rho to each clinic based on intervention arm
         Ps = ifelse(TRTArm == 1, ratep[2], ratep[1])) # assign base rate of outcome to each clinic based on intervention arm

# compare Ps to rateP2; observed versus true

probclinic <- (rnorm(n_site, 0, sqrt(sigma_site))) # create random effects for clinic level 
prob_hcl <- rnorm(ngroup1, 0, sqrt(sigma_hcl)) # create random effects for HCS level

Data$ClinicEff <- c(probclinic); Data$SystemEff = c(A2%*%prob_hcl) # assigning the random normal variates to each clinic and HCS

Data <- Data %>% mutate(Ps = c(expit(logit(Ps) + ClinicEff + SystemEff))) # rewrite over Ps for clinic specific probability of vaccination
# evaluate the probability of vaccination for each clinic

Data <- Data %>% rowwise() %>% mutate(Y_Obs = rbinom(1,ceiling(Frequency), Ps)) # draw observed vaccinated count

# calculate true Y count
Data$Y_true = apply(mapply(function(x,y,z) rmultinom(1,size = x, prob = c(z, (y*(1-z)), (1 - z - y*(1-z)))),
                           Data$Frequency, Data$Rhos, Data$Ps), 2, function(q) q[1] + q[2]) 

Y = Data$Y_Obs; N_x = ceiling(Data$Frequency)
Y_true = Data$Y_true; N_True = Data$N_true

###########################################################################################################

# Fitting GLMER to observed data 

mod1 <- lme4::glmer(cbind(Y,N_x - Y) ~ X[,2] + (1|id2)+(1|id1), family= 'binomial')
summary(mod1)

# polya-gamma with fuzzy outcomes data

#nsamp = 2000; warmup = 500 # number of posterior samples drawn.

pg1 <-  site_pgamma4(Data_Frames = Data_Frames,Y=Y,N=N_x,nsamps=2000,warmup=500, hyp = c(0.01, 0.01), underreport = TRUE,
                           Beta_Params=Beta_Set[[1]], 
                           Csi_Params = csi1_and2, Csi_Params2 = csi1_and2)
# test out the model with the first set of specifications (Beta_Set[[1]])
# you can change the 1 to 2 or 3 to analyze the data set with the two other specifications

betasampsd <- apply(pg1$Beta_Samps[[1]],1,function(x) expit(sum(x)))
betasampsd0 <- expit((pg1$Beta_Samps[[1]])[,1])
res3 <- c(mean(betasampsd), quantile(betasampsd, probs = c(0.025, 0.975)), 
          mean(betasampsd0), quantile(betasampsd0, probs = c(0.025, 0.975)), 
          colMeans(pg1$Beta_Samps[[1]])[2], quantile(pg1$Beta_Samps[[1]][,2], probs = c(.025, 0.975)),
          mean(exp(pg1$Beta_Samps[[1]][,2])), quantile(exp(pg1$Beta_Samps[[1]][,2]), probs = c(.025, 0.975)),
          median(1/pg1$Tau_samps[[1]]),median(1/pg1$Tau_samps[[2]]))

displayres <- data.frame("Estimates" = res3[c(1,4,7,10)], "Lower 95% CI" = res3[c(2, 5, 8,11)], "Upper 95% CI" = res3[c(3,6,9,12)],
                         "Parameters" = c("Outcome Rate in Intervention Arm", "Outcome Rate in Control Arm",
                                          "Log-Odds Ratio of Intervention", "Odds Ratio of Intervention"))
varestims <- res3[c(13,14)]
names(varestims) <- c("Median Estimate for Clinic Level Variance", "Median Estimate for HCS Level Variance")
# in current analysis: we are under-specifying the misclassified outcomes in the treatment arm 
# therefore the odds ratio we estimate is higher than the observed odds ratio but still less than truth
displayres; varestims
