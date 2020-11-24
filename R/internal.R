#' @useDynLib gfiExtremes
#' @importFrom Rcpp evalCpp
#' @import coda
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

#' @importFrom stats optim
#' @noRd
selectThreshold <- function(X, candidates){
  params <- matrix(NA_real_, nrow = length(candidates), ncol = 2L)
  values <- rep(NA_real_, length(candidates))
  for(i in seq_along(candidates)){
    mu <- candidates[i]
    fn <- function(gamma_sigma){
      - sum(dgpareto(
        x = X[X >= mu], 
        mu = mu, 
        gamma = gamma_sigma[1L]/(1-gamma_sigma[1L]), 
        sigma = gamma_sigma[2L]/(1-gamma_sigma[2L]), 
        log = TRUE
      ))
    }
    opt <- optim(
      par = c(0.5, 0.5), fn = fn, 
      method = "L-BFGS-B", lower = 0.01, upper = 0.99, 
      control = list(maxit = 500L)
    )
    if(opt[["convergence"]] == 0L){
      params[i, ] <- opt[["par"]]/(1-opt[["par"]])
      values[i, ] <- opt[["value"]]
    }
  }
  imin <- which.min(values)
  c(mu = candidates[imin], gamma = params[imin, 1L], sigma = params[imin, 2L])
}