# ripw
An R package for Reshaped Inverse Propensity Weighted (RIPW) Two-way fixed effects (TWFE) estimator (authored by Xiaoman Luo and Lihua Lei)

## Overview
This R package implements the RIPW-TWFE estimator for panel data analysis proposed in our paper: [Design-Robust Two-Way-Fixed-Effects Regression For Panel Data](https://arxiv.org/abs/2107.13737). 

- `ripw_twfe()` produces the RIPW-TWFE estimator for any doubly average treatment effect (DATE) with any user-specified weights, reshaped distribution, treatment model estimator, and outcome model estimator, with and without derandomized cross-fitting. 
- `muhat_twfe()` produces $(\hat{\mu}_{it}(1), \hat{\mu}_{it}(0))$ by fitting the TWFE estimator with covariate-treatment interactions described in Section 3.2. In particular, it allows for time-varying/invariant covariates with time-varying/invariant coefficients. It is called when no user-defined outcome model estimator is given for `ripw_twfe()`.
- `pihat_cox()` produces the estimated generalized propensity score $\hat{\pi}_i(w)$ by fitting the Cox model described in Section 5.2. In particular, it allows for time-varying/invariant covariates with time-varying/invariant coefficients. It is called when no user-defined treatment model estimator is given for `ripw_twfe()` and the treatment is staggered. 
- `solve_date_maxmin` produces the solution of a DATE equation, if any, that maximizes $\min_{w}\Pi(w)$. The procedure is described in Appendix C.5. It is called when no user-defined reshaped distribution is given for `ripw_twfe()`.
- `data(opentable)` loads into R the OpenTable dataset analyzed in Section 5.2. 

## Installation

```
if (!require("devtools")){
    install.packages("devtools")
}
devtools::install_github("lihualei71/ripw")
```

We suggest installing [survival](https://cran.r-project.org/web/packages/survival/index.html) and [lfe](https://cran.r-project.org/web/packages/lfe/index.html) to take advantage of the built-in learners. To run the examples, [tidyverse](https://www.tidyverse.org/packages/) needs to be installed. 

## Usage Examples
We illustrate the usage of ripw package on a simulated panel data with staggered treatment.
```
# Load package
library("ripw")

# Generate data
set.seed(1)
n <- 1000
T <- 4
d_tv <- 2 # number of time-varying covariates
d_ti <- 5 # number of time-invariant covariates

Y <- matrix(rnorm(n*T), nrow = n)
adopt_time <- sample(0:T, n, replace = TRUE) # generate adoption times
tr <- sapply(adopt_time, function(x){
	if (x == 0){
		rep(0, T)
	} else {
		c(rep(0, x-1), rep(1,T-x+1))
	}
})
tr <- t(tr)
X_tv <- lapply(1:T, function(x){
	mat <- matrix(rnorm(n*d_tv), nrow = n)
	colnames(mat) <- paste0("Xtv", 1:d_tv)
	mat
}) # generate time-varying covariates as a list of nxd_tv matrices
X_ti <- matrix(rnorm(n*d_ti), nrow = n) # generate time-invariant covariates as an nxd_ti matrix
colnames(X_ti) <- paste0("Xti", 1:d_ti)

# Set up the parameters needed for treatment and outcome modeling
# set of variables for pihat_cox with time-varying coefficients
pi_tvc <- c("Xtv1", "Xti1") 
# set of variables for pihat_cox with time-invariant coefficients
pi_tic <- c("Xtv2", "Xti2", "Xti3", "Xti4") 
pihat_params <- list(X_tv = X_tv, X_ti = X_ti, 
                     tvc = pi_tvc, tic = pi_tic)

# set of variables for muhat_twfe with time-varying coefficients 
# in the non-interacted terms
mu_main_tvc <- c("Xtv1", "Xti1") 
# set of variables for muhat_twfe with time-invariant coefficients
# in the non-interacted terms; no time-invariant variable should 
# be included since they will be collinear with unit fixed effects
mu_main_tic <- c("Xtv2") 
# set of variables for muhat_twfe with time-varying coefficients 
# in the interacted terms
mu_int_tvc <- c("Xtv1", "Xti1", "Xti2") 
# set of variables for muhat_twfe with time-invariant coefficients
# in the interacted terms
mu_int_tic <- c("Xtv1", "Xtv2", "Xti1", "Xti3", "Xti4") 
# joint=TRUE means fit a joint TWFE regression on both treated and control units
muhat_params <- list(X_tv = X_tv, X_ti = X_ti, 
                     main_tvc = mu_main_tvc, main_tic = mu_main_tic, 
                     int_tvc = mu_int_tvc, int_tic = mu_int_tic, 
					 joint = TRUE) 

# Fit the RIPW-TWFE estimator
nreps <- 20 # number of sample splits for derandomization
nfolds <- 10 # number of folds for cross-fitting
res <- ripw_twfe(Y, tr, 
                 pihat_params = pihat_params, 
				 muhat_params = muhat_params, 
				 nreps = nreps, nfolds = nfolds)
```


