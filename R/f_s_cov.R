#' f_* covariance function
#'
#' Calculates the posterior covariance matrix for test input X_* for a
#' given draw from the posterior of w.
#'
#' @param w_p A draw from the posterior of w
#' @param k_xx K(X,X) matrix (from the training data)
#' @param k_sx K(X*, X) matrix
#' @param k_xs K(X, X*) matrix
#' @param k_ss K(X*,X*) matrix
#' @param x Matrix of training covariates
#'
#' @export
f_s_cov <- function(w_p, k_ss, k_sx, k_xs, k_xx, x) {
    chol_matr <- chol(k_xx + diag(as.vector(exp(t(w_p) %*% t(x)))))
    f_s_cov <- k_ss - k_sx %*% solve(chol_matr) %*% solve(t(chol_matr)) %*% k_xs
    return(f_s_cov)
}