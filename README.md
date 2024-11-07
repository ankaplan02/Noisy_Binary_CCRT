# Group-Level Bayesian Logistic Regression with Negative Misclassification of Binary Outcomes

This addition to Bayesian logistic regression attends to when binary outcomes and their corresponding sizes (i.e., $`Y,N`$) are misrepresented by the available data in the context of group-level analysis (e.g., cluster randomized trials). Let $`Y_{s}`$ and $`N_{s}`$ denote the observed counts of success and size for group $s$. The error we assume is that we are undercounting the successes and overcounting the number eligible. 

The key example situation is vaccinations where the analysis is taken over a study period. The primary outcome may be any vaccination during the study period however the eligibility to be in the study is that the patients cannot be vaccinated. With imperfect record keeping, we may have strong evidence to suggest that we have missed true study-period related vaccinations while we overlooked eligibility (or patients didn't report they were previously vaccinated). 

All but three steps in this method reflect that of a hierarchical generalized linear mixed model for binary outcomes, which is commonly used in cluster randomized trials. The additional steps are the following theoretical corrections to the data: 

- $`N_{s}^{*} = N_{s} - \lceil \rho_{s,2}(N_{s} - Y_{s}\rceil - \lceil \rho_{s,3}Y_{s}\rceil`$ (excluding those erroneously labeled as eligible in the study)
- $`E_{s} = Y_{s} - \lceil \rho_{s,3}Y_{s}\rceil`$ (eligible observed success count)
- $`Y_{s}^{*} = E_{s} + \lceil \rho_{s,1}(N_{s}^{*} - E_{s})\rceil`$ (reclassifying eligible unvaccinated)

Prior elicitation is required to fulfill these steps. A prior distribution for the vector $`\rho_{s} = [\rho_{s,1}, \rho_{s,2}, \rho_{s,3}]`$ needs to be specified for each site $s$. The most complexity we expect is that misclassification rates differ by intervention arm; however, direct specification of each site's misclassification rate vector can be done with an example in the available simulation code. We specify that for a generic $\rho$, that it follows a truncated Beta distribution. This is how we solve for its specifications. 

We assume the prior distribution for $\rho$, denoted as $`\pi(\rho)`$, is equal to $`\text{Beta}(\gamma, \lambda)_{[a,b]}`$, also known as the truncated Beta distribution. The values of $a$ and $b$ are the lower and upper limits we can specify for plausible values for the misclassification rate $\rho$.

The research team can easily specify $a$ and $b$, but the values of $\gamma$ and $\lambda$ require more information. We can estimate values for these two parameters provided a most likely value for $\rho$ (i.e., the mode), and lower- and upper-end percentiles for $\rho$, denote these as $\kappa_{p_{1}}$ and $\kappa_{p_{2}}$, respectively. Under $\pi(\rho)$, the mode of $\rho$, $M= \frac{\gamma-1}{\gamma+\lambda-2}$; rearranging yields $\lambda = \frac{\gamma(1-M) + (2M -1)}{M}$. Last, using $\kappa_{p_{1}}$ and $\kappa_{p_{2}}$ we can further refine the shape of the desired Beta distribution by minimizing the following objective function: 

```math
\Bigg[\int_{\kappa_{p_{2}}}^{b}\text{Beta}(\gamma, \lambda)_{[a,b]}d\rho -
\int_{a}^{\kappa_{p_{1}}}\text{Beta}(\gamma, \lambda)_{[a,b]}d\rho\Bigg] - (p_{2} - p_{1})=0.
```
The important features needed to optimize this function are the following: 

- $M$ the mode, the most likely value we think $\rho$ will take on. 
- $a$ and $b$, the lower and upper limits we will allow $\rho$ to plausibly take on
- and finally the shape defining parameters, $\kappa_{p_{1}}$ and $\kappa_{p_{2}}$, where in the code will need to be separated like so: c($`\kappa_{p_{1}}, p_{1}, \kappa_{p_{2}}, p_{2}`$); for analogy, one can specify that $p_{1}=0.025$ and $p_{2} = 0.975$ for a "95% confidence-like interval" for the value of $\rho$, where the limits of this interval are $\kappa_{0.025}$ and $\kappa_{0.975}$. This was done in the associated methods article.

Once these values are specified and input into the function in the correct way, the algorithm will find the appropriate $\gamma$ and $\lambda$ values. Ensure that $a < \kappa_{p_{1}} < \kappa_{p_{2}} < b$ and $p_{1} < p_{2}$ such that a solution to the above can be attained. 

Once these prior distribution assumptions are incorporated, the rest follows Bayesian Logistic Regression with random effects but the data is perterbed with each posterior sample by the distributions we set for $\rho_{s}$. 
