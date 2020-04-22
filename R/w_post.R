#' Posterior distribution for w
#'
#' A slightly longer description here.
#'
#' @param w Values of the paramaters
#' @param dist_args Vector of function arguments, in order they are
#' y, x, length_scale, sigma_f, and tau2.
#' Written as list(y, x, length_scale, sigma_f, tau2).
#'
#' @export
w_post <- function(w, dist_args) {

    y <- dist_args[[1]]
    x <- dist_args[[2]]
    length_scale <- dist_args[[3]]
    sigma_f <- dist_args[[4]]
    tau2 <- dist_args[[5]]

    k <- se_kern(x_1 = x, x_2 = x, length_scale, sigma_f)
    l <- chol(k + diag(as.vector(exp(t(w) %*% t(x)))))
    log_det <- 2 * sum(log(diag(l)))

    post_density <- -0.5 * log_det -
        0.5 * t(y) %*% solve(l) %*% solve(t(l)) %*% y -
        t(w) %*% w / (2 * tau2)
    return(post_density)
}