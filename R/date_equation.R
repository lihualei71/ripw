#### Non-linear solver for the DATE equation (Appendix C.5)

#' Solver for the DATE equation with \eqn{\min_{w}\Pi(w) \ge \Pi_{low}} (Appendix C.5)
#'
#' @param support: the matrix \eqn{(\tilde{w}_{(1)}, ..., \tilde{w}_{(K)})}
#' @param xi: weights used for the DATE
#' @param init: initial value
#' @param Pi_low: lower bound for min_{w}Pi(w)
#' @return
#' \item{sol}{a solution (if found)}
#' \item{val}{the value of the objective at sol}
#' \item{init_val}{the value of the objective at init}
#'
#' @noRd
solve_date <- function(support, xi,
                       init = rep(1, ncol(support)) / ncol(support),
                       Pi_low = 0){
    T <- nrow(support)
    K <- ncol(support)
    A <- support
    JA <- A - rep(1, T) %*% t(colMeans(A))
    bmat <- A * JA - xi %*% t(diag(t(JA) %*% A))
    dateeq <- function(p){
        Ap <- A %*% p
        LHS <- bmat %*% p
        RHS <- (diag(as.numeric(Ap)) - xi %*% t(Ap)) %*% (Ap - mean(Ap))
        LHS - RHS
    }
    init_val <- dateeq(init)
    if (max(abs(init_val)) < 1e-12){
        cat("Initial value is a solution\n")
        list(sol = obj$par, val = init_val, init_val = init_val)
    }
    normal_fac <- max(1e-6, sum(init_val^2) / 2)
    fn <- function(p){
        sum(dateeq(p)^2) / 2 / normal_fac
    }
    gr <- function(p){
        Ap <- A %*% p
        JAp <- Ap - mean(Ap)
        comp1 <- (diag(as.numeric(Ap)) - xi %*% t(Ap)) %*% JA
        comp2 <- (diag(as.numeric(JAp)) - xi %*% t(JAp)) %*% A
        Dp <- bmat - comp1 - comp2
        t(Dp) %*% dateeq(p) / normal_fac
    }
    ui1 <- rep(1, K)
    ui2 <- rep(-1, K)
    ui3 <- diag(rep(1, K))
    ui <- rbind(ui1, ui2, ui3)
    ci <- c(1-1e-4, -1-1e-4, rep(max(Pi_low - 1e-4, 0), K))

    obj <- constrOptim(init, fn, gr, ui, ci, outer.eps = 1e-6)
    val <- dateeq(obj$par)
    list(sol = obj$par, val = val, init_val = init_val)
}

#' Solver for the DATE equation that maximizes min_{w}Pi(w) (Appendix C.5)
#'
#' @param support the matrix \eqn{(\tilde{w}_{(1)}, ..., \tilde{w}_{(K)})}
#' @param xi weights used for the DATE
#' @param init initial value
#' @param tol tolerance for the objective value
#'
#' @return
#' \item{sol}{a solution (if found)}
#' \item{val}{the value of the objective at sol}
#' \item{init_val}{the value of the objective at init}
#'
#' @export
solve_date_maxmin <- function(support, xi,
                              init = rep(1, ncol(support)) / ncol(support),
                              tol = 1e-6){
    fn <- function(Pi_low){
        obj <- solve_date(support, xi, init, Pi_low)
        max(abs(obj$val))
    }
    if (fn(0) > tol){
        cat("No solution is found\n")
        return(NULL)
    }
    left <- 0
    right <- min(init)
    if (fn(right) < tol){
        return(solve_date(support, xi, init, right))
    }
    while (right - left > 1e-4){
        mid <- (left + right) / 2
        if (fn(mid) > tol){
            right <- mid
        } else {
            left <- mid
        }
    }
    return(solve_date(support, xi, init, left))
}

#' Find the type of treatment paths
#'
#' @noRd
treatment_type <- function(tr){
    T <- ncol(tr)
    support <- unique(tr)
    if (all(support[, -1] - support[, -T] >= 0)){
        return("staggered")
    } else if (all(rowSums(support) <= 1)){
        return("transient")
    } else {
        return("other")
    }
}

#' Default reshaped distribution
#'
#' @noRd
default_Pi_fun <- function(tr, xi){
    T <- ncol(tr)
    type <- treatment_type(tr)
    if (type == "staggered" && sum(abs(xi - rep(1/T, T))) < 1e-6) {
        Pi_fun <- function(w){
            if (sum(w) %in% c(0, T)){
                (T + 1) / 4 / T
            } else {
                1 / 2 / T
            }
        }
    } else if (type == "transient" && sum(abs(xi - rep(1/T, T))) < 1e-6){
        Pi_fun <- function(w){
            1 / (T+1)
        }
    } else {
        warning("No default treatment type detected. Solving the DATE equation using the generic BFGS algorithm to maximize min_w Pi(w).")
        support <- t(unique(tr))
        Pi <- solve_date_maxmin(support, xi)$sol
        Pi_fun <- function(w){
            id <- which.min(colMeans((support - w)^2))
            Pi[id]
        }
        # stop("No default reshaped distribution is found")
    }
    return(Pi_fun)
}

#' Default pihat function
#'
#' @noRd
default_pihat_fun <- function(tr){
    type <- treatment_type(tr)
    if (type == "staggered"){
        return(pihat_cox)
    } else {
        stop("No default pihat function is found")
    }
}
