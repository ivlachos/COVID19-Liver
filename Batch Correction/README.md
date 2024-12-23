Building on approaches that use residuals from a negative binomial generalized linear model (NB-GLM) to normalize single cell data, 
we fitted a NB-GLM using an efficient implementation of a Gamma-Poisson GLM with batch as the covariates. 
We then used the deviance residuals from this model as the expression adjusted for batch effects. 
For downstream analysis that required counts, we also generated counts corrected for batch by expanding and scaling the model described by Zhang et al. 
using a scalable implementation of a Gamma-Poisson GLM.
