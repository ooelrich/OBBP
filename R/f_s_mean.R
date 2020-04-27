#' f_* mean function
#'
#' Calculates the posterior mean for test input X_* for a given draw
#' from the posterior of w.
#'
#' @param w_p A draw from the posterior of w
#' @param k_xx K(X,X) matrix (from the training data)
#' @param k_sx K(X*, X) matrix
#' @param y Observed outcomes
#' @param x Matrix of training covariates
#'
#' @export
f_s_mean <- function(w_p, k_xx, k_sx, y, x) {
    chol_matr <- chol(k_xx + diag(as.vector(exp(t(w_p) %*% t(x)))))
    f_s_mean <- k_sx %*% solve(chol_matr) %*% solve(t(chol_matr)) %*% y
    return(f_s_mean)
}