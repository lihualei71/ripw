check_staggered <- function(tr){
    row_staggered <- apply(tr, 1, function(w){
        all(diff(w) >= 0)
    })
    all(row_staggered)
}

cox_standardize_input <- function(X_tv = NULL, X_ti = NULL,
                                  T = NULL,
                                  tvc = NULL, tic = NULL){
    ## Return an empty matrix if nothing is given
    if ((is.null(X_tv) && is.null(X_ti)) ||
        (is.null(tvc) && is.null(tic))){
        X <- NULL
        formula <- "survival::Surv(time1, time2, status)~1"
        return(list(X = NULL, formula = formula,
                    vars = NULL))
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
    X_tvc <- matrix(NA, nrow = nrow(X[[1]])*T, ncol = 0)
    if (!is.null(tvc)){
        X_tvc <- lapply(X, function(x){
            x[, tvc, drop=FALSE]
        })
        X_tvc <- do.call(Matrix::bdiag, X_tvc)
        X_tvc <- as.matrix(X_tvc)
        tvc <- paste0(rep(tvc, T),
                           ".T",
                           rep(1:T, each = length(tvc)))
        colnames(X_tvc) <- tvc
    }

    ## Construct the matrix for main effects with time-invariant coeffcients
    X_tic <- matrix(NA, nrow = nrow(X[[1]])*T, ncol = 0)
    if (!is.null(tic)){
        X_tic <- lapply(X, function(x){
            x[, tic, drop=FALSE]
        })
        X_tic <- do.call(rbind, X_tic)
        colnames(X_tic) <- tic
    }

    ## Combine two covariate matrices
    X <- cbind(X_tvc, X_tic)

    ## Generate variable names for main effects and interactions
    vars <- c(tvc, tic)

    ## Generate formula
    formula <- paste0("survival::Surv(time1, time2, status)~",
                      paste0(vars,
                             collapse = "+"))

    return(list(X = X, formula = formula,
                vars = vars))
}

#' Cross-fit a Cox model to estimate the generalized propensity score for staggered treatments
#'
#' \code{pihat_cox} fits a Cox model with time-varying/invariant variables and time-varying/invariant coefficients.
#' To formally write down the model, we define the following quantities
#' \itemize{
#' \item{\eqn{X_{itj}}}{ j-th time-varying covariate \eqn{(j=1,\ldots,p)} for unit i at time t \eqn{(j=p+1,\ldots, p+q)}}
#' \item{\eqn{X_{ik}}}{ j-th time-invariant covariate for unit i \eqn{(j=p+1,\ldots, p+q)}}
#' \item{\eqn{S_{tvc}\subset \{1, \ldots, p+q\}}}{ the subset of variables with time-varying coefficients}
#' \item{\eqn{S_{tic}\subset \{1, \ldots, p+q\}}}{ the subset of variables with time-invariant coefficients}
#' }
#' Then the Cox model we fit is
#' \deqn{\lambda(t\mid X_i) = \lambda_0(t)\exp\Bigg\{\sum_{j\in S_{tvc}, j\le p}X_{itj}\beta_{jt} + \sum_{j\in S_{tvc}, j > p}X_{ij}\beta_{jt}
#' + \sum_{j\in S_{tic}, j\le p}X_{itj}\beta_{j} + \sum_{j\in S_{tic}, j > p}X_{ij}\beta_{j}\Bigg\}.}
#'
#'
#' @param tr an nxT binary matrix for the treatment variable (which must be staggered)
#' @param X_tv a length-T list of nxp matrices where p is the number of time-varying variables. The t-th element is the matrix \eqn{(X_{itj})_{i\in [n], j\in [p]}}
#' @param X_ti an nxq matrix where q is the number of time-invariant variables
#' @param tvc a vector of variable names for (both time-varying and time-invariant) variables with time-varying coefficients. Leave it as default if no time-varying coefficient is used.
#' @param tic a vector of variable names for (both time-varying and time-invariant) variables with time-invariant coefficients. Leave it as default if no time-invariant coefficient is used.
#' @param cf TRUE iff cross-fitting is used
#' @param nfolds number of folds (K in Algorithm 1)
#' @param foldid NULL if fold indices are generated randomly or a list of length K with the k-th element a vector of indices in fold k
#'
#' @return a length-n vector of generalized propensity score estimates
#'
#' @examples
#' \donttest{
#' # Load opentable data
#' data(opentable)
#' library("tidyverse")
#'
#' # Generate the treatment assignments
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
#' tic <- c("confirmed", colnames(X_ti))
#' gps <- pihat_cox(tr, X_tv, X_ti, tic = tic, cf = TRUE, nfolds = nfolds)
#'
#' # Using cross-fitting and randomly generated folds
#' set.seed(2023)
#' foldid <- ripw:::gen_cf_inds(nrow(tr), nfolds)
#' gps2 <- pihat_cox(tr, X_tv, X_ti, tic = tic, foldid = foldid)
#' identical(gps, gps2)
#'
#' # Not using cross-fitting
#' gps3 <- pihat_cox(tr, X_tv, X_ti, tic = tic, cf = FALSE)
#' }
#'
#' @export
pihat_cox <- function(tr,
                      X_tv = NULL, X_ti = NULL,
                      tvc = NULL, tic = NULL,
                      cf = TRUE, nfolds = 10,
                      foldid = NULL){
    ## Generate data
    n <- nrow(tr)
    T <- ncol(tr)
    data <- data.frame(tr = as.numeric(tr))
    data$unit <- rep(1:n, T)
    data$time <- rep(1:T, each = n)

    ## Check is the assignments are staggered
    if (!check_staggered(tr)){
        stop("Assignments are not staggered")
    }

    ## Check foldid
    if (!is.null(foldid)){
        nfolds <- length(foldid)
    } else {
        foldid <- gen_cf_inds(n, nfolds)
    }

    ## Standardize X
    obj <- cox_standardize_input(X_tv, X_ti, T,
                                 tvc, tic)

    ## Combine X with tr
    data <- data.frame(data, obj$X)
    if (ncol(data) > 3){
        names(data)[4:ncol(data)] <- colnames(obj$X)
    }

    ## Generate data for coxph
    survdata <- list()
    adopt_time <- rep(0, n) # 0 for never-treated units
    unit <- NULL # to pass the R CMD check
    for (i in 1:n){
        tmpdata <- data[data$unit == i, ]
        inds <- which(tmpdata$tr > 0)
        if (length(inds) > 0){
            ind <- min(inds)
            tmpdata <- tmpdata[1:ind, ]
            tmpdata$status <- c(rep(0, ind-1), 1)
            adopt_time[i] <- ind
        } else {
            tmpdata$status <- 0
        }
        tmpdata$time1 <- tmpdata$time-1
        tmpdata$time2 <- tmpdata$time
        survdata[[i]] <- tmpdata
    }
    survdata <- do.call(rbind, survdata)

    ## Fit cox model
    gps_est <- rep(NA, n)
    for (k in 1:nfolds){
        ids <- foldid[[k]]
        train <- survdata[!(survdata$unit %in% ids), ]
        test <- survdata[survdata$unit %in% ids, ]
        fit <- survival::coxph(as.formula(obj$formula), data = train)
        res <- summary(survival::survfit(fit, newdata = test, id = unit))
        for (index in ids){
            adopt <- adopt_time[index] > 0 # if ever-treated
            if (is.null(res$strata)){
                survprob <- res$surv
            } else {
                survprob <- res$surv[res$strata == index]
            }
            if (length(survprob) == 0){
                ## If adopt_time[index] is earlier than
                ## adoption times of all training units,
                ## use P(T <= start_date) instead
                earliest_date <- min(summary(survival::survfit(fit))$time)
                newtest <- test[test$unit == index, ]
                added_row <- newtest[nrow(newtest), ]
                added_row$time2 <- earliest_date
                added_row$time1 <- earliest_date - 1
                added_row$status <- 1
                newtest <- rbind(newtest, added_row)
                newres <- summary(survival::survfit(fit, newdata = newtest, id = unit))
                survprob <- newres$surv
            }
            if (adopt){
                if (length(survprob) == 1){
                    gps_est[index] <- 1 - survprob
                } else {
                    gps_est[index] <- -diff(tail(survprob, 2))
                }
            } else {
                ## For never-treated units, we regularize
                ## the gps as P(T >= end_date - 1)
                gps_est[index] <- tail(survprob, 2)[1]
            }
        }
    }
    return(gps_est)
}
