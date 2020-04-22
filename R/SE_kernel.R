SE_kernel <-
function(X, lengthScale, sigma_f){
    K <- matrix(NA, nrow = nrow(X), ncol = nrow(X))
    for(i in 1:nrow(X)){
        for(j in 1:nrow(X)){
            K[i, j] <- sigma_f * exp(-sum((X[i, 2] - X[j, 2])^2) / (2*lengthScale))
        }
    }
    return(K)
}
