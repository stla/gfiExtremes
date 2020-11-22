#' Generalized Pareto distribution
#' @description Density, distribution function, quantile function, and random 
#'   generation for the generalized Pareto distribution.
#'   
#' @param x numeric vector
#' @param q numeric vector of quantiles
#' @param p numeric vector of probabilities
#' @param n positive integer, the desired number of simulations
#' @param mu location parameter
#' @param gamma shape parameter, positive
#' @param sigma scale parameter, strictly positive
#' @param log logical, whether to return the log-density
#' 
#' @importFrom stats runif
#' 
#' @rdname GPareto
#' @name GPareto
#' @export
dgpareto <- function(x, mu, gamma, sigma, log = FALSE){ 
  stopifnot(gamma >= 0, sigma > 0)
  out <- if(log) log(numeric(length(x))) else numeric(length(x))
  less_than_mu <- q < mu
  in_support <- !less_than_mu
  if(any(in_support)){
    x <- x[in_support]
    z <- (x - mu) / sigma
    out[in_support] <- if(gamma == 0){
      if(log) -z else exp(-z)
    }else{
      if(log){
        (-1/gamma-1) * log1p(gamma*z) - log(sigma)
      }else{
        (1 + gamma*z)^(-1/gamma-1) / sigma
      }
    }
  }
  out
}

#' @rdname GPareto
#' @export
pgpareto <- function(q, mu, gamma, sigma){ 
  stopifnot(gamma >= 0, sigma > 0)
  out <- numeric(length(q))
  less_than_mu <- q < mu
  in_support <- !less_than_mu
  if(any(in_support)){
    q <- q[in_support]
    z <- (q - mu) / sigma
    out[in_support] <- if(gamma == 0){
      1 - exp(-z)
    }else{
      1 - (1 + gamma*z)^(-1/gamma)
    }
  }
  out
}

#' @rdname GPareto
#' @export
rgpareto <- function(n, mu, gamma, sigma){ 
  qgpareto(runif(n), mu, gamma, sigma)
}

#' @rdname GPareto
#' @export
qgpareto <- function(p, mu, gamma, sigma){
  stopifnot(all(p >= 0 & p <= 1))
  stopifnot(gamma >= 0, sigma > 0)
  if(gamma == 0){
    mu - sigma * log1p(-p)
  }else{
    mu + sigma * ((1-p)^(-gamma) - 1) / gamma
  }
}
