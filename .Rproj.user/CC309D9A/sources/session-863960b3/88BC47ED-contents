#### RIPW-TWFE estimator (Algorithm 1)

#' RIPW-TWFE estimator and inference with derandomized cross-fitting
#'
#' \code{ripw_twfe} fits the two-way-fixed-effects regression with modified outcomes weighted by the ratio
#' between the reshaped distribution and the generalized propensity score. It allows for user-defined estimators of
#' treatment model (i.e., \eqn{\hat{\pi}_i}) and outcome model (i.e., \eqn{\hat{\mu}_{it}(1)} and \eqn{\hat{\mu}_{it}(0)})
#'
#' \code{Pi_fun} is a function that maps a treatment path to its density under the reshaped distribution. \code{Pi_fun}
#' should take a single input as a binary vector that gives the treatment path and return a non-negative scalar that gives
#' the density. The default value of \code{Pi_fun} is NULL, in which case the internal function \code{default_Pi_fun} will
#' be used. \code{default_Pi_fun} first checks if the treatment is staggered, transient (one-off), or other. For staggered
#' treatment, if \code{xi = rep(1/T, T)}, it uses the solution (C.5) in the paper, i.e., the never-treated and
#' always-treated get (T+1)/4T and other assignments get 1/2T. For transient treatment, if \code{xi = rep(1/T, T)}, it uses
#' the uniform distribution over the never-treated and all ever-treated assignments. In all other cases, it calls
#' \code{\link{solve_date_maxmin}} and wraps the output into a function.
#'
#' \code{pihat_fun} is a function that generates the estimates of generalized propensity scores. It must output a nx1 vector
#' of generalized propensity scores for all units. The inputs must include \code{cf} and \code{foldid} (see
#' \code{\link{pihat_cox}} for details) and can include other arguments passed through \code{pihat_params}.
#'
#' \code{muhat_fun} is a function that generates the estimates of\eqn{\hat{\mu}_{it}(1)} and \eqn{\hat{\mu}_{it}(0)}. It must
#' output a \code{list(muhat1 = , muhat0 = )} where each element is a nxT matrix. The inputs must include \code{cf}
#' and \code{foldid} (see \code{\link{muhat_twfe}} for details) and can include other arguments passed through \code{muhat_params}.
#'
#'
#' @param Y an nxT matrix for the outcome variable
#' @param tr an nxT binary matrix for the treatment variable
#' @param xi DATE weights; a T-dim vector of nonnegative numbers that sum up to 1
#' @param Pi_fun a function to generate the reshaped distribution; see Details
#' @param pihat_fun a function to estimate generalized propensity scores; see Details
#' @param muhat_fun a function to estimate regression adjustments; see Details
#' @param pihat_params a list of extra parameters for pihat_fun
#' @param muhat_params a list of extra parameters for muhat_fun
#' @param cf TRUE iff cross-fitting is used
#' @param nfolds number of folds (K in Algorithm 1)
#' @param foldid NULL if fold indices are generated randomly or a list of length K with the k-th element a vector of indices in fold k
#' @param nreps number of data splits (B in Algorithm 1)
#' @param Theta_bound a truncation bound for \eqn{\Theta_i \triangleq \Pi(w) / \hat{\pi}_i(w)}; \code{Inf} by default
#' @param show_progress_bar TRUE if progress bar is shown. FALSE by default if \code{nreps=1}.
#'
#' @return
#' \item{tauhat}{point estimate}
#' \item{se}{standard error}
#' \item{tstat}{t statistic}
#' \item{pval}{p-value for the two-sided test of zero effect}
#'
#' @examples
#' \donttest{
#' # Load opentable data
#' data(opentable)
#' library("tidyverse")
#'
#' # Generate the outcomes and treatment assignments
#' Y <- opentable %>% select(state, time, reserv.diff) %>%
#'  spread(time, reserv.diff) %>%
#'  select(-state) %>%
#'  as.matrix
#' tr <- opentable %>% select(state, time, treat) %>%
#'  spread(time, treat) %>%
#'  select(-state) %>%
#'  as.matrix
#' T <- ncol(Y)
#'
#' ## Construct covariate matrices
#' # vote, beds, and region are time-invariant covariates
#' X_vote <- opentable$vote %>%
#'   matrix(nrow = T) %>% t %>% .[, 1]
#' X_beds <- opentable$beds %>%
#'   matrix(nrow = T) %>% t %>% .[, 1]
#' X_region <- opentable$region %>%
#'   matrix(nrow = T) %>% t %>% .[, 1] %>%
#'   data.frame(region = .) %>%
#'   model.matrix(~0+region, data = .) %>%
#'   .[, -1] %>% as.matrix
#' X_ti <- cbind(vote = X_vote, beds = X_beds, X_region)
#'
#' # confirmed is a time-varying covariate
#' X_confirmed <- opentable$confirmed %>%
#'   matrix(nrow = T) %>% t
#' X_tv <- lapply(1:T, function(t){
#'   tmp <- as.matrix(X_confirmed[, t])
#'   colnames(tmp) <- "confirmed"
#'   tmp
#' })
#'
#' # Derandomized cross-fitting
#' set.seed(2023)
#' nfolds <- 10
#' nreps <- 10
#'
#' ## We will fit a Cox model for generalized propensity scores and a TWFE regression with interactions
#' ## for the regression adjustment.
#' # X_tv is for time-varying covariates, X_ti is for time-invariant covariates
#' # tic is the list of variables with time-invariant coefficients; see pihat_cox for details
#' pihat_params <- list(X_tv = X_tv,
#'                     X_ti = X_ti,
#'                     tic = c("confirmed", colnames(X_ti)))
#' # main_tic is the list of variables with time-invariant coefficients in the main regression,
#' # int_tic is the list of variables with stime-invariant coefficients in the interaction terms;
#' # see muhat_twfe for details
#' muhat_params <- list(X_tv = X_tv,
#'                      X_ti = X_ti,
#'                      main_tic = "confirmed",
#'                      int_tic = c("confirmed", colnames(X_ti)),
#'                      joint = TRUE)
#'
#' # Compute the derandomized cross-fitted RIPW-TWFE estimator with K=10 and B=10
#' res <- ripw_twfe(Y, tr,
#'                  pihat_fun = pihat_cox,
#'                  muhat_fun = muhat_twfe,
#'                  pihat_params = pihat_params,
#'                  muhat_params = muhat_params,
#'                  nreps = nreps,
#'                  nfolds = nfolds)
#' print(res)
#' }
#'
#' @export
ripw_twfe <- function(Y, tr,
                      xi = rep(1/ncol(Y), ncol(Y)),
                      Pi_fun = NULL,
                      pihat_fun = NULL,
                      muhat_fun = NULL,
                      pihat_params = list(),
                      muhat_params = list(),
                      cf = TRUE, nfolds = 10,
                      foldid = NULL,
                      nreps = 1,
                      Theta_bound = Inf,
                      show_progress_bar = ifelse(nreps == 1, FALSE, TRUE)){
    n <- nrow(Y)
    if (is.null(Pi_fun)){
        Pi_fun <- default_Pi_fun(tr, xi)
    }
    Pi <- apply(tr, 1, Pi_fun)

    if (is.null(pihat_fun)){
        pihat_fun <- default_pihat_fun(tr)
    }

    if (is.null(muhat_fun)){
        muhat_fun <- muhat_twfe
    }

    if (!is.null(foldid)){
        nreps <- 1
    } else {
        foldid <- gen_cf_inds(n, nfolds)
    }

    pihat_params$tr <- tr
    muhat_params$Y <- Y
    muhat_params$tr <- tr

    denom <- rep(NA, nreps)
    tauhat <- rep(NA, nreps)
    vhat <- matrix(NA, n, nreps)

    if (show_progress_bar){
        pb <- txtProgressBar(min = 0, max = nreps, style = 3)
    }

    for (i in 1:nreps){
        ## Estimate pi
        pihat <- do.call(pihat_fun,
                         c(pihat_params, list(cf = TRUE, nfolds = 10, foldid = foldid)))
        Theta <- Pi / pihat
        Theta <- pmin(Theta, Theta_bound)
        Theta <- Theta / mean(Theta) # normalize Theta

        ## Estimate mu
        muhat <- do.call(muhat_fun,
                         c(muhat_params, list(cf = TRUE, nfolds = 10, foldid = foldid)))
        nuhat <- muhat$muhat1 - muhat$muhat0
        mhat <- muhat$muhat0
        for (ids in foldid){
            mhat[ids, ] <- center(mhat[ids, ])
            nuhat[ids, ] <- nuhat[ids, ] - sum(colMeans(nuhat[ids, ]) * xi)
        }

        ## Generate modified outcomes
        Ytd <- Y - mhat - nuhat * tr

        ## Center Ytd and tr to get JY and JW
        Ytd_c <- Ytd - rowMeans(Ytd)
        tr_c <- tr - rowMeans(tr)

        Gamma_w <- colMeans(Theta * tr_c)
        Gamma_y <- colMeans(Theta * Ytd_c)
        Gamma_ww <- mean(Theta * rowSums(tr_c^2))
        Gamma_wy <- mean(Theta * rowSums(tr_c * Ytd_c))
        denom[i] <- Gamma_ww - sum(Gamma_w^2)
        numer <- Gamma_wy - sum(Gamma_w * Gamma_y)
        tauhat[i] <- numer / denom[i]

        Ytd_c <- Ytd_c - tr_c * tauhat[i]
        vhat[, i] <- (Gamma_wy - tauhat[i] * Gamma_ww + rowSums(tr_c * Ytd_c) - Ytd_c %*% Gamma_w - tr_c %*% (Gamma_y - Gamma_w * tauhat[i])) * Theta

        ## Generate folds for next iteration
        foldid <- gen_cf_inds(n, nfolds)

        if (show_progress_bar){
            setTxtProgressBar(pb,i)
        }
    }

    close(pb)
    tauhat <- sum(tauhat * denom) / sum(denom)
    v <- rowSums(vhat) / sum(denom)
    tau_se <- sd(v) / sqrt(n)
    tstat <- tauhat / tau_se
    pval <- pnorm(abs(tstat), lower.tail = FALSE) * 2

    return(list(tauhat = tauhat, se = tau_se,
                tstat = tstat, pval = pval))
}
