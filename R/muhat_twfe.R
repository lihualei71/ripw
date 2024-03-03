twfe_standardize_input <- function(X_tv = NULL, X_ti = NULL,
                                   T = NULL,
                                   main_tvc = NULL, main_tic = NULL,
                                   int_tvc = NULL, int_tic = NULL,
                                   joint = TRUE){
    ## No interaction is allowed for separate fitting
    if (!joint && (!is.null(int_tvc) || !is.null(int_tic))){
        warning("No interaction is included when joint=FALSE")
    }

    ## Combine time-varying and time-invariant covariates
    if (is.null(X_tv)){
        X <- lapply(1:T, function(t){
            X_ti
        })
    } else if (!is.null(X_ti)){
        X <- lapply(X_tv, function(x){
            cbind(x, X_ti)
        })
    } else {
        X <- X_tv
    }
    T <- length(X)

    ## Construct the matrix for main effects with time-varying coeffcients
    X_main_tvc <- matrix(NA, nrow = nrow(X[[1]])*T, ncol = 0)
    if (!is.null(main_tvc)){
        X_main_tvc <- lapply(X, function(x){
            x[, main_tvc, drop=FALSE]
        })
        X_main_tvc <- do.call(Matrix::bdiag, X_main_tvc)
        X_main_tvc <- as.matrix(X_main_tvc)
        main_tvc <- paste0(rep(main_tvc, T),
                           ".T",
                           rep(1:T, each = length(main_tvc)))
        colnames(X_main_tvc) <- main_tvc
    }

    ## Construct the matrix for main effects with time-invariant coeffcients
    X_main_tic <- matrix(NA, nrow = nrow(X[[1]])*T, ncol = 0)
    if (!is.null(main_tic)){
        X_main_tic <- lapply(X, function(x){
            x[, main_tic, drop=FALSE]
        })
        X_main_tic <- do.call(rbind, X_main_tic)
        colnames(X_main_tic) <- main_tic
    }

    ## Construct the matrix for interactive effects with time-varying coeffcients
    X_int_tvc <- matrix(NA, nrow = nrow(X[[1]])*T, ncol = 0)
    if (!is.null(int_tvc)){
        X_int_tvc <- lapply(X, function(x){
            x[, int_tvc, drop=FALSE]
        })
        X_int_tvc <- do.call(Matrix::bdiag, X_int_tvc)
        X_int_tvc <- as.matrix(X_int_tvc)
        int_tvc <- paste0(rep(int_tvc, T),
                          ".T",
                          rep(1:T, each = length(int_tvc)),
                          "_tr")
        colnames(X_int_tvc) <- int_tvc
    }

    ## Construct the matrix for interactive effects with time-invariant coeffcients
    X_int_tic <- matrix(NA, nrow = nrow(X[[1]])*T, ncol = 0)
    if (!is.null(int_tic)){
        X_int_tic <- lapply(X, function(x){
            x[, int_tic, drop=FALSE]
        })
        X_int_tic <- do.call(rbind, X_int_tic)
        int_tic <- paste0(int_tic, "_tr")
        colnames(X_int_tic) <- int_tic
    }

    ## Combine four covariate matrices
    if (joint){
        X <- cbind(X_main_tvc, X_main_tic, X_int_tvc, X_int_tic)
    } else {
        X <- cbind(X_main_tvc, X_main_tic)
    }

    ## Generate formula
    if (joint){
        formula <- paste0("y~",
                          paste0(c(main_tvc, main_tic,
                                   int_tvc, int_tic),
                                 collapse = "+"),
                          "+tr|unit+time")
    } else {
        formula <- paste0("y~",
                          paste0(c(main_tvc, main_tic),
                                 collapse = "+"),
                          "|unit+time")
    }

    ## Generate variable names for main effects and interactions
    main <- c(main_tvc, main_tic)
    interact <- c(int_tvc, int_tic)

    return(list(X = X, formula = formula,
                main = main, interact = interact))
}

## Estimate mhat from a two-way fixed effect regression
## without cross-fitting
##
## Add unit fixed effects estimates into \hat{m}_{0} if usefe=TRUE
## when no cross-fitting is used. Does not make a difference for
## RIPW but does so for aggregated AIPW
twfe_muhat_onefold <- function(data,
                               formula,
                               Xtest,
                               main, interact, joint){
    options(warn = -1)

    ## Generate interaction terms
    if (!is.null(interact) && joint){
        data[, interact] <- data[, interact] * data$tr
    }

    ## Fit the regression models
    if (!joint){
        fit1 <- lfe::felm(formula, data = data, subset = (data$tr == 1))
        fit0 <- lfe::felm(formula, data = data, subset = (data$tr == 0))
        coef1 <- as.numeric(coef(fit1))
        coef1[is.nan(coef1)] <- 0
        coef0 <- as.numeric(coef(fit0))
        coef0[is.nan(coef0)] <- 0
        muhat1 <- Xtest %*% coef1
        muhat0 <- Xtest %*% coef0
    } else {
        fit <- lfe::felm(formula, data = data)
        coef <- coef(fit)
        coef_main <- coef[main]
        coef_main[is.nan(coef_main)] <- 0
        muhat0 <- Xtest[, main, drop=FALSE] %*% as.matrix(coef_main)
        if (!is.null(interact) ){
            ## The variable name could either be "tr:xxx" or "xxx:tr"
            coef_tau <- coef[interact]
            coef_tau[is.nan(coef_tau)] <- 0
            tauhat <- Xtest[, interact, drop=FALSE] %*% coef_tau
            muhat1 <- muhat0 + tauhat
        } else {
            muhat1 <- muhat0
        }
    }
    muhat1 <- matrix(muhat1, ncol = T)
    muhat0 <- matrix(muhat0, ncol = T)

    options(warn = 0)
    muhat <- list(muhat1 = muhat1, muhat0 = muhat0)
    return(muhat)
}

#' Estimate \eqn{\hat{\mu}_{it}(1)} and \eqn{\hat{\mu}_{it}(0)} from a TWFE regression with cross-fitting
#'
#' \code{muhat_twfe} fits a TWFE regression time-varying/invariant variables and time-varying/invariant coefficients in both non-interacted and interacted terms.
#' To formally write down the model, we define the following quantities
#' \itemize{
#' \item{\eqn{X_{itj}}}{ j-th time-varying covariate \eqn{(j=1,\ldots,p)} for unit i at time t \eqn{(j=p+1,\ldots, p+q)}}
#' \item{\eqn{X_{ik}}}{ j-th time-invariant covariate for unit i \eqn{(j=p+1,\ldots, p+q)}}
#' \item{\eqn{S_{main\_tvc}\subset \{1, \ldots, p+q\}}}{ the subset of variables with time-varying coefficients in non-interacted terms}
#' \item{\eqn{S_{main\_tic}\subset \{1, \ldots, p+q\}}}{ the subset of variables with time-invariant coefficients in non-interacted terms}
#' \item{\eqn{S_{int\_tvc}\subset \{1, \ldots, p+q\}}}{ the subset of variables with time-varying coefficients in interacted terms}
#' \item{\eqn{S_{int\_tic}\subset \{1, \ldots, p+q\}}}{ the subset of variables with time-invariant coefficients in interacted terms}
#' }
#' When \code{joint=TRUE}, we fit the following TWFE regression on the treated and control groups together:
#' \deqn{Y_{it} = \alpha_i + \lambda_t + \sum_{j\in S_{main\_tvc}, j\le p}X_{itj}\beta_{jt} + \sum_{j\in S_{main\_tvc}, j > p}X_{ij}\beta_{jt}
#' + \sum_{j\in S_{main\_tic}, j\le p}X_{itj}\beta_{j} + \sum_{j\in S_{main\_tic}, j > p}X_{ij}\beta_{j} + \tau W_{it} + \sum_{j\in S_{int\_tvc}, j\le p}W_{it}X_{itj}\phi_{jt}
#' + \sum_{j\in S_{int\_tvc}, j > p}W_{it}X_{ij}\phi_{jt} + \sum_{j\in S_{int\_tic}, j\le p}W_{it}X_{itj}\phi_{j} + \sum_{j\in S_{int\_tic}, j > p}W_{it}X_{ij}\phi_{j} + \epsilon_{it}.}
#' With this fit,
#' \deqn{\hat{\mu}_{it}(0) =  \sum_{j\in S_{main\_tvc}, j\le p}X_{itj}\beta_{jt} + \sum_{j\in S_{main\_tvc}, j > p}X_{ij}\beta_{jt}
#' + \sum_{j\in S_{main\_tic}, j\le p}X_{itj}\beta_{j} + \sum_{j\in S_{main\_tic}, j > p}X_{ij}\beta_{j},}
#' and
#' \deqn{\hat{\mu}_{it}(1) = \hat{\mu}_{it}(0) + \sum_{j\in S_{int\_tvc}, j\le p}X_{itj}\phi_{jt}
#' + \sum_{j\in S_{int\_tvc}, j > p}X_{ij}\phi_{jt} + \sum_{j\in S_{int\_tic}, j\le p}X_{itj}\phi_{j} + \sum_{j\in S_{int\_tic}, j > p}X_{ij}\phi_{j}.}
#' When \code{joint=FALSE}, we fit the following TWFE regression for the treated and control groups separately:
#' \deqn{Y_{it} = \alpha_i + \lambda_t + \sum_{j\in S_{main\_tvc}, j\le p}X_{itj}\beta_{jt} + \sum_{j\in S_{main\_tvc}, j > p}X_{ij}\beta_{jt}
#' + \sum_{j\in S_{main\_tic}, j\le p}X_{itj}\beta_{j} + \sum_{j\in S_{main\_tic}, j > p}X_{ij}\beta_{j} + \epsilon_{it}.}
#' In particular, \code{main_tic} and \code{int_inc} will not be used.
#‘
#' @param Y an nxT matrix for the outcome variable
#' @param tr an nxT binary matrix for the treatment variable (which must be staggered)
#' @param X_tv a length-T list of nxp matrices where p is the number of time-varying variables. The t-th element is the matrix \eqn{(X_{itj})_{i\in [n], j\in [p]}}
#' @param X_ti an nxq matrix where q is the number of time-invariant variables
#' @param main_tvc a vector of variable names for (both time-varying and time-invariant) variables with time-varying coefficients in the non-interacted terms. Leave it as default if no time-varying coefficient is used.
#' @param main_tic a vector of variable names for (both time-varying and time-invariant) variables with time-invariant coefficients in the non-interacted terms. Leave it as default if no time-invariant coefficient is used.
#' @param int_tvc a vector of variable names for (both time-varying and time-invariant) variables with time-varying coefficients in the interacted terms. Leave it as default if no time-varying coefficient is used.
#' @param int_tic a vector of variable names for (both time-varying and time-invariant) variables with time-invariant coefficients in the interacted terms. Leave it as default if no time-invariant coefficient is used.
#' @param joint TRUE if treatment and control groups are jointly fit. See Details.
#' @param cf TRUE iff cross-fitting is used
#' @param nfolds number of folds (K in Algorithm 1)
#' @param foldid NULL if fold indices are generated randomly or a list of length K with the k-th element a vector of indices in fold k
#'
#' @return a list of form list(muhat1 = , muhat0 = ) where muhat1 and muhat0 are both nxT matrices
#'
#' @examples
#' \donttest{
#' # Load opentable data
#' data(opentable)
#' library("tidyverse")
#'
#' # Generate the treatment assignments
#' Y <- opentable %>% select(state, time, reserv.diff) %>%
#'  spread(time, reserv.diff) %>%
#'  select(-state) %>%
#'  as.matrix
#' tr <- opentable %>% select(state, time, treat) %>%
#'  spread(time, treat) %>%
#'  select(-state) %>%
#'  as.matrix
#' T <- ncol(tr)
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
#' # Using cross-fitting and randomly generated folds
#' set.seed(2023)
#' nfolds <- 10
#' main_tic <- "confirmed"
#' int_tic <- c("confirmed", colnames(X_ti))
#' muhat <- muhat_twfe(Y, tr, X_tv, X_ti, main_tic = main_tic, int_tic = int_tic,
#'                     joint = TRUE, cf = TRUE, nfolds = nfolds)
#'
#' # Using cross-fitting and randomly generated folds
#' set.seed(2023)
#' foldid <- ripw:::gen_cf_inds(nrow(Y), nfolds)
#' muhat2 <- muhat_twfe(Y, tr, X_tv, X_ti, main_tic = main_tic, int_tic = int_tic,
#'                      joint = TRUE, foldid = foldid)
#' identical(muhat, muhat2)
#'
#' # Not using cross-fitting
#' muhat3 <- muhat_twfe(Y, tr, X_tv, X_ti, main_tic = main_tic, int_tic = int_tic,
#'                      joint = TRUE, cf = FALSE)
#' }
#' @export
muhat_twfe <- function(Y, tr,
                       X_tv = NULL,
                       X_ti = NULL,
                       main_tvc = NULL, main_tic = NULL,
                       int_tvc = NULL, int_tic = NULL,
                       joint = TRUE,
                       cf = TRUE, nfolds = 10,
                       foldid = NULL){
    ## Check the format of Y and tr
    n <- nrow(Y)
    T <- ncol(Y)
    if (nrow(tr) != n || ncol(tr) != T){
        stop("dim(Y) does not match dim(tr)")
    }

    ## Return zeros when no covariate is included
    if ((is.null(X_tv) && is.null(X_ti)) ||
        (is.null(main_tvc) && is.null(main_tic) && is.null(int_tvc) && is.null(int_tic))){
        muhat1 <- muhat0 <- matrix(0, n, T)
        return(list(muhat1 = muhat1, muhat0 = muhat0))
    }

    ## Transform Y, and tr into nT X 1 vectors
    Y <- as.numeric(Y)
    tr <- as.numeric(tr)

    ## Add columns for unit and time FE
    data <- data.frame(y = Y, tr = tr)
    data$unit <- as.factor(rep(1:n, T))
    data$time <- as.factor(rep(1:T, each = n))

    ## Standardize X
    obj <- twfe_standardize_input(X_tv, X_ti, T,
                                  main_tvc, main_tic,
                                  int_tvc, int_tic,
                                  joint)

    ## Combine X with Y and tr
    data <- data.frame(data, obj$X)
    if (ncol(data) > 4){
        names(data)[5:ncol(data)] <- colnames(obj$X)
    }

    ## Check foldid
    if (!is.null(foldid)){
        nfolds <- length(foldid)
    } else {
        foldid <- gen_cf_inds(n, nfolds)
    }

    ## Run twfe_muhat_onefold if nfolds=1
    if (!cf || nfolds == 1){
        muhat <- twfe_muhat_onefold(
            data, as.formula(obj$formula), obj$X,
            obj$main, obj$interact, joint)
        muhat$muhat1 <- matrix(muhat$muhat1, ncol = T)
        muhat$muhat0 <- matrix(muhat$muhat0, ncol = T)
        return(muhat)
    }

    ## Cross-fitting
    muhat1 <- muhat0 <- matrix(NA, n, T)
    for (i in 1:nfolds){
        inds <- foldid[[i]]
        data_inds <- which(data$unit %in% inds)
        data_train <- data[-data_inds, ,drop=FALSE]
        Xtest <- obj$X[data_inds, ,drop=FALSE]
        muhat <- twfe_muhat_onefold(
            data_train, as.formula(obj$formula), Xtest,
            obj$main, obj$interact, joint)
        muhat1[inds, ] <- muhat$muhat1
        muhat0[inds, ] <- muhat$muhat0
    }
    muhat <- list(muhat1 = muhat1, muhat0 = muhat0)
    return(muhat)
}

