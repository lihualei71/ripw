## Generate indices for cross-fitting
gen_cf_inds <- function(n, nfolds){
    tmp_inds <- sample(1:n, n)
    foldid <- quantile(1:n, probs = 0:nfolds / nfolds, type = 1)
    foldid[1] <- 0
    foldid <- lapply(1:nfolds, function(i){
        sort(tmp_inds[(foldid[i]+1):foldid[i+1]])
    })
    return(foldid)
}

## Center rows and columns of vectors and matrices
center <- function(dat){
    if (is.vector(dat) || (is.matrix(dat) && ncol(dat) == 1)){
        dat <- dat - mean(dat)
    } else if (is.matrix(dat)){
        rowavg <- rowMeans(dat)
        dat <- dat - rowavg
        colavg <- colMeans(dat)
        dat <- t(t(dat) - colavg)
    } else {
        stop("dat must be a vector or a matrix")
    }
    return(dat)
}
