# This algorithm was proposed in the paper 'Robust Computation of Linear Models by Convex
# Relaxation' GM function returns a geometric median.

GM <- function(X) {
    delta <- 1e-10
    epsilon <- 1e-10
    Loop <- 100
    
    data <- as.matrix(X)
    N <- nrow(data)
    d <- ncol(data)
    cn <- rep(0, d)
    cs <- rep(0, d)
    cn <- as.vector(apply(data, 2, function(v) mean(v)))
    for (l in 1:Loop) {
        ds <- sapply(1:N, function(i) 1/max(delta, norm(data[i, ] - cn, type = "2")))
        cs <- (t(data) %*% ds)/sum(ds)
        
        if (norm(cs - cn, type = "2") < epsilon) {
            cn <- cs
            break
        } else {
            cn <- cs
        }
    }
    return(cn)
}
