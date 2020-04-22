#' metropolis_sampler
#'
#' This function implements a simple Metropolis MC function
#'
#' @param startVals Starting values for the parameters
#' @param y Response variable
#' @param X matrix of covariates
#' @param iter Number of iterations the sampler will run
#' @param burnin Number of iterations that will be discarded
#'
#' @export

metropolis_sampler <- function(startVals, y, X, iter, burnin) {
    reject <- 0
    theta <- matrix(NA, ncol = 2, nrow = iter)
    theta[1, ] <- startVals
    for (i in 2:iter) {
        prop <- mvrnorm(1, theta[i - 1, ],
                        matrix(c(0.3, -0.1, -0.1, 0.3),
                        nrow = 2))
        if(postdist(prop, y, X) - postdist(theta[i - 1, ], y, X) > log(runif(1))) {
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
