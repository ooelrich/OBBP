#' Metropolis MC sampler
#'
#' This function implements a simple Metropolis MCMC algorithm.
#'
#' @param dist Distribution to draw from
#' @param dist_args List of arguments needed for the posterior, in this
#' case it should be list(y = responde, x = covariates, length_scale, sigma_f, tau2)
#' @param start_vals Starting values for the parameters
#' @param iter Number of iterations the sampler will run
#' @param burnin Number of iterations that will be discarded
#'
#' @export
metropolis_sampler <- function(dist, dist_args, start_vals, iter, burnin) {

    reject <- 0
    theta <- matrix(NA, ncol = length(start_vals), nrow = iter)
    theta[1, ] <- start_vals

    for (i in 2:iter) {
        prop <- MASS::mvrnorm(1, theta[i - 1, ],
                        matrix(c(0.3, -0.1, -0.1, 0.3),
                        nrow = 2))
        lr <- dist(prop, dist_args) - dist(theta[i - 1, ], dist_args)
        if (lr > log(stats::runif(1))) {
            theta[i, ] <- prop
        }
        else{
            theta[i, ] <- theta[i - 1, ]
            if (i > burnin) {
                reject <- reject + 1
            }
        }

    }

    print(reject)
    print(iter - burnin)
    return(theta[(burnin + 1):iter, ])
}