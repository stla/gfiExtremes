#' Title
#'
#' @param X numeric vector of data
#' @param beta vector of probabilities corresponding to the quantiles to be 
#'   estimated
#' @param threshold.init a guess of the unknown threshold, must be in the range 
#'   of \code{X}
#' @param gamma.init starting value for gamma in the MCMC
#' @param sigma.init starting value for sigma in the MCMC
#' @param sd.gamma standard deviation for the proposed gamma
#' @param sd.sigma standard deviation for the proposed sigma
#' @param p1 probability that the MCMC will propose a new \code{(gamma,sigma)}; 
#'   \code{(1-p1)} would be the probability that the MCMC chain will propose a 
#'   new index for a new threshold
#' @param p2 probability that the new index proposed will be larger than the 
#'   current index
#' @param lambda1 the small jump the index variable will make
#' @param lambda2 the large jump the index variable will make; happens 1 of 
#'   every 10 iterations
#' @param Jnumb number of subsamples that are taken from the Jacobian
#' @param iter number of iterations per chain (burnin excluded)
#' @param burnin number of the first MCMC iterations discarded
#' @param thin thinning number for the MCMC chain. (e.g. if it is 1 no iteration 
#'   is skipped)
#' @param nchains number of MCMC chains to run
#' @param nthreads number of threads to run the chains in parallel
#' @param seeds the seeds used for the MCMC sampler; one seed per chain, or 
#'   \code{NULL} to use random seeds
#'
#' @return An object of class \code{\link[coda:mcmc]{mcmc}} if \code{nchains=1}, 
#'   otherwise an object of class \code{\link[coda:mcmc.list]{mcmc.list}}.
#'   
#' @references Damian V. Wandler & Jan Hannig. 
#'   \emph{Generalized fiducial confidence intervals for extremes}.
#'   Extremes (2012) 15:67–87.
#'   <doi:10.1007/s10687-011-0127-9>
#' 
#' @export
#' @importFrom ismev gpd.fit
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach `%dopar%`
#'
#' @examples set.seed(666L)
#' X <- rgamma(500L, shape = 10, rate = 1)
#' gf <- gfigpd2(X, beta = c(0.98, 0.99), iter = 3000L)
#' summary(gf)
#' qgamma(c(0.98, 0.99), shape = 10, rate = 1)
#' thresholdEstimate(gf)
#' rejectionRate(gf)
#' HPDinterval(gf)
#' HPDinterval(joinMCMCchains(gf))
gfigpd2 <- function(
  X, beta, threshold.init = NA, 
  gamma.init = NA, sigma.init = NA, sd.gamma = NA, sd.sigma = NA, 
  p1 = 0.9, p2 = 0.5, lambda1 = 2, lambda2 = 10, Jnumb = 50L, 
  iter = 10000L, burnin = 2000L, thin = 6L,
  nchains = nthreads, nthreads = parallel::detectCores(), seeds = NULL) {
  
  stopifnot(thin >= 1L, nchains >= 1L, nthreads >= 1L)
  nthreads <- min(nthreads, nchains)
  
  X <- sort(X) # -->> so there's no need to sort in C++

  if(is.na(threshold.init)){
    i <- floor(0.85 * length(X))
    threshold.init <- X[i]
  }else{
    if(threshold.init <= X[1L] || threshold.init >= X[length(X)]){
      stop(
        "The value of `threshold.init` is not in the range of `X`."
      )
    }
    i <- thresholdIndex(X, threshold.init)
  }
  
  # Initialize the default values for the tuning parameters of the MCMC chain
  if(is.na(gamma.init) || is.na(sigma.init)) {
    mle.fit <- gpd.fit(X, X[i], show = FALSE)
    if(is.na(gamma.init)) gamma.init <- mle.fit$mle[2L]
    if(is.na(sigma.init)) sigma.init <- mle.fit$mle[1L]
  }
  if(is.na(sd.gamma)) sd.gamma <- 2 * abs(gamma.init) / 3
  if(is.na(sd.sigma)) sd.sigma <- 2 * sigma.init / 3
  
  skip.number <- thin - 1L
  number.iterations <- (skip.number + 1L) * iter + burnin
  
  if(is.null(seeds)){
    seed1 <- sample.int(2000000L, 1)
    seeds <- seed1 + 2000000L * (0L:(nchains-1L))
  }else{
    if(length(seeds) != nchains){
      stop(
        "Please specify one seed per chain."
      )
    }
    seeds <- abs(as.integer(seeds))
  }
  
  params <- c("gamma", "sigma", "index", paste0("beta", seq_along(beta)))
  
  # run the MCMC chain
  if(nchains == 1L){
    chain <- thinChain(MCMCchain(
      X, beta, gamma.init, sigma.init, threshold.init, i, 
      p1, p2, lambda1, lambda2, sd.gamma, sd.sigma,
      number.iterations, burnin, Jnumb, seeds[1L]
    ), skip.number)
    colnames(chain) <- params
    chain[, 4L:(3L + length(beta))] <- chain[, 4L:(3L + length(beta))] + X[1L]
    thresholdEstimate <- X[median(chain[,"index"])]
  }else{
    if(nthreads == 1L){
      chains <- vector("list", nchains)
      for(k in 1L:nchains){
        chains[[k]] <- MCMCchain(
          X, beta, gamma.init, sigma.init, i, 
          p1, p2, lambda1, lambda2, sd.gamma, sd.sigma,
          number.iterations, burnin, Jnumb, seeds[k]
        )
      }
    }else{
      # nblocks <- ceiling(nchains/nthreads)
      # blocks <- vector("list", nblocks)
      # nthreads <- 
      #   c(rep(nthreads, nblocks-1L), nchains - ((nblocks-1L)*nthreads))
      # for(b in 1L:nblocks){
      #   registerDoParallel(cores = nthreads[b])
      # }
      cl <- makeCluster(nthreads)
      registerDoParallel(cl)
      chains <- foreach(
        k = 1L:nchains, .combine = list, .multicombine = TRUE, 
        .export = "MCMCchain"
      ) %dopar% MCMCchain(
        X, beta, gamma.init, sigma.init, threshold.init, i, 
        p1, p2, lambda1, lambda2, sd.gamma, sd.sigma,
        number.iterations, burnin, Jnumb, seeds[k]
      )
      stopCluster(cl)
    }
    chains <- lapply(chains, thinChain, skip = skip.number)
    chains <- lapply(chains, `colnames<-`, value = params)
    chains <- lapply(chains, function(chain){
      chain[, 4L:(3L + length(beta))] <- chain[, 4L:(3L + length(beta))] + X[1L]
      chain
    })
    index <- c(vapply(chains, function(chain){
      chain[, "index"]
    }, FUN.VALUE = numeric(iter)))
    thresholdEstimate <- X[median(index)]
  }
  
  if(nchains == 1L){
    out <- coda::mcmc(chain, start = burnin+1L, thin = thin)
  }else{
    out <- 
      coda::mcmc.list(lapply(chains, coda::mcmc, start = burnin+1L, thin = thin))
  }
  
  attr(out, "beta") <- beta
  attr(out, "threshold") <- thresholdEstimate
  
  return(out)

  # # Indictor for the acceptance rate. // pas besoin avec coda::rejectionRate
  # ii <- integer(nrow(x.t))
  # for (i in 1L:(nrow(x.t) - 1L)) {
  #   if (x.t[i, 4L] != x.t[i + 1L, 4L]) {
  #     ii[i + 1L] <- 1L
  #   }
  # }
  # acceptance.rate <- mean(ii)
  # cat("acceptance rate: ", acceptance.rate)
  # 
  # x.t[, 4L:(3L + length(beta))] <- x.t[, 4L:(3L + length(beta))] + X[1L]
  # 
  # return(x.t)
}
