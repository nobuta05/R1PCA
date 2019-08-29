curr_dir <- getwd()
frame_files <- lapply(sys.frames(), function(x) x$ofile)
frame_files <- Filter(Negate(is.null), frame_files)
setwd(dirname(frame_files[[length(frame_files)]]))
source("GeometricMedian.r")
setwd(curr_dir)

r1pca <- function(x, maxd, centered=1, isCentered=FALSE) {
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

  inpX <- as.matrix(x)
  N <- nrow(inpX)
  d <- ncol(inpX)
  
  # centered==1の時は`\min_{c\in\mathbb{R}^d} \sum_{i=1}^{N} \left\| x_i - c \right\|_{2}` を最小にする`\hat{c}\in\mathbb{R}^d` を中心ベクトル(平均ベクトル)として中心化を行う
  # centered==2の時は通常の平均ベクトルによって中心化を行う
  mu <- rep(0, d)
  if(!isCentered && centered == 1) {
    mu <- GM(inpX)
  }
  else if(!isCentered && centered == 2) {
    mu <- as.vector(apply(data, 2, function(v) mean(v)))
  }
  X <- inpX - matrix(mu, N, d, T)

  U <- (eigen(t(X) %*% X)$vectors)[,1:maxd]
  s <- sapply(
    1:N,
    function(i) sqrt( max(0, ( X[i,] %*% ((diag(d) - U%*% t(U))%*%X[i,]) )[1]) )
  )
  c <- median(s)

  for(l in 1:Loop) {
    q <- sapply(
      1:N,
      function(i) wH(X[i,], U, c)
    )
    W <- t(X) %*% (X * q)
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

  scores <- X %*% U
  return(list(center=mu, loadings=U, scores=scores))
}
