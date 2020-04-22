#' Squared exponential kernel
#'
#' An implementation of the squared exponential kernel, K(x_1, x_2).
#' Note that x_1 = x_2 is valid.
#'
#' @param x_1 First design matrix
#' @param x_2 Second design matrix
#' @param length_scale The length scale
#' @param sigma_f Some kind of noise
#'
#' @export
se_kern <- function(x_1, x_2, length_scale, sigma_f) {

    if (ncol(x_1) != ncol(x_2)) {
        stop("x_1 and x_2 must have the same number of columns")
    }

    kern <- matrix(NA, nrow = nrow(x_1), ncol = nrow(x_2))

    for (i in seq_len(nrow(x_1))) {
        for (j in seq_len(nrow(x_2))) {
            kern[i, j] <- sigma_f * exp(-sum((x_1[i, ] - x_2[j, 2])^2) / 
                            (2 * length_scale))
        }
    }

    return(kern)
}