#' Posterior distrbution
#' 
#' @param w values of the paramaters
#' @param y Response variable
#' @param X Covariates
#' 
#' @export

postdist <- function(w, y, X){
    L <- chol(K + diag(as.vector(exp(t(w) %*% t(X)))))
    logDet <- 2*sum(log(diag(L)))
    likelihood <- -0.5 * logDet - 
        0.5 * t(y) %*% solve(L) %*% solve(t(L)) %*% y - 
        t(w) %*% w / (2 * tau2) 
    return(likelihood)
}
