source("GeometricMedian.R")

r1pca <- function(x, k, centered=2) {
  epsilon <- 1e-10
  delta <- 1e-10
  Loop <- 100
  wH <- function(v, U, c) {
    if(norm(v-U%*%t(U)%*%v, type="2") <= c) {
      return(1)
    }
    else {
      return(c/norm(v-U%*%t(U)%*%v, type="2"))
    }
  }

  data <- as.matrix(x)
  N <- nrow(data)
  d <- ncol(data)
  if (k>d)
    k <- d

  # centered==1の時は`\min_{c\in\mathbb{R}^d} \sum_{i=1}^{N} \left\| x_i - c \right\|_{2}` を最小にする`\hat{c}\in\mathbb{R}^d` を中心ベクトル(平均ベクトル)として中心化を行う
  # centered==2の時は通常の平均ベクトルによって中心化を行う
  mu <- rep(0, d)
  if(centered == 1) {
    mu <- GM(data)
    data <- data - matrix(t(mu), N, d, T)
  }
  else if(centered == 2) {
    mu <- as.vector(apply(data, 2, function(v) mean(v)))
    data <- data - matrix(t(mu), N, d, T)
  }

  U <- (eigen(t(data) %*% data)$vectors)[,1:k]
  s <- sapply(
    1:N,
    function(i) sqrt( max(0, ( data[i,] %*% ((diag(d) - U%*% t(U))%*%data[i,]) )[1]) )
  )
  c <- median(s)

  for(l in 1:Loop) {
    q <- sapply(
      1:N,
      function(i) wH(data[i,], U, c)
    )
    W <- t(data) %*% (data * q)
    WUqr <- qr(W %*% U)
    Us <- qr.Q(WUqr)
    if(norm(Us - U) < epsilon) {
      U <- Us
      break
    }
    else {
      U <- Us
    }
  }

  scores <- data %*% U
  return(
    list(
      center=mu,
      loadings=U,
      scores=scores
    )
  )
}
