#' Newton's method for finding optima
#'
#' A simple implementation of Newton's method to find posterior mode.
#' For use with the Laplace approximation in exercise 2, lab1.
#'
#' @param y Vector of outcome variables.
#' @param k_xx Covariance matrix of the GP.
#' @param precision Measured by the squared distance between two consequtive
#' f vectors. Defaults to 0.05.
#'
#' @export
newtons_method <- function(y, k_xx, precision = 0.05) {

    f <- rep(0, length(y))
    f_old <- rep(1, length(y))
    
    while ((t(f-f_old) %*% (f-f_old)) > precision) {
        w <- diag(exp(as.vector(f)))
        ch <- chol(solve(k_xx + diag(rep(0.01, 16))) + w)
        f_old <- f
        f <- solve(ch) %*% solve(t(ch)) %*% (w %*% f + y - exp(f))
    }

    return(f)
}
