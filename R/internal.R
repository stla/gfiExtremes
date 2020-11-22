#' @useDynLib gfiExtremes
#' @importFrom Rcpp evalCpp
NULL

rGPDmodel <- function(nsims, mu = 10, gamma = 1, sigma = 1){
  U <- runif(nsims)
  P <- qgpareto(1-U, mu, gamma, sigma)
  p <- 1 / (1/mu + 1)   # p * 1/a = 1 - p => p = 1/(1/a+1)
  X <- numeric(nsims)
  for(i in 1L:nsims){
    X[i] <- sample(c(mu*U[i], P[i]), 1L, prob = c(p, 1-p))
  }
  X
}

thinChain <- function(chain, skip){
  niterations <- nrow(chain)
  every.ith <- c(TRUE, rep(FALSE, skip))
  keep <- rep(every.ith, ceiling(niterations / (skip+1L)))[1L:niterations]
  chain[keep, ]
}

thresholdIndex <- function(Xs, a){
  # Xs is X sorted
  match(TRUE, Xs - a >= 0)
}
