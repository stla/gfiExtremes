library(ismev)

######## This is the actual function to run for the package############


Fiducial.GPD.unknown.thresh <- function(
  X, beta, CI.level, g = NaN, s = NaN, i = NaN, p1 = .9, p2 = .5, 
  lambda1 = 2, lambda2 = 10, sd.g = NaN, sd.s = NaN, skip.number = NaN, 
  chain.length = 10000, burnin = 2000, J.numb = 50) {
  # X is the vector of data.
  # beta is the value or vector of to estimate the Beta-quantile.
  # CI.level is the confidence level for confidence intervals.
  # g is the starting value for gamma in the MCMC chain.
  # s is the starting value for sigma in the MCMC chain.
  # i is the indix for the intial threshold at the X(i) order statistic.
  # p1 is the probability that the MCMC will propose a new (gamma,sigma).
  # (1-p1) would be the probability that the MCMC chain will propose a new index for a new threshold.
  # p2 is the probability that the new index proposed will be larger than the current index.
  # lambda1 is the small jump the index variable will make.
  # lambda2 is the large jump the index variable will make. Happens 1 of every 10 iterations.
  # sd.g is the standard deviation for the proposed gamma.
  # sd.s is the standard deviation for the proposed sigma.
  # skip.number is a thinning number for the MCMC chain. (e.g. if skip.number is 1 we will use ever other MCMC iteration).
  # chain.length is the length of the MCMC chain.
  # burnin is number of the first MCMC iterations omitted from the chain.
  # J.numb is the number of subsamples that are taken from the Jacobian.
  
  X <- sort(X)
  n <- length(X)
  
  # Intialize the default values for the tuning parameters of the MCMC chain.
  
  if (i == "NaN") i <- n * .85
  
  mle.fit <- gpd.fit(X, X[i], show = FALSE)
  if (g == "NaN") g <- mle.fit$mle[2]
  if (s == "NaN") s <- mle.fit$mle[1]
  
  if (sd.g == "NaN") sd.g <- 2 * abs(g) / 3
  if (sd.s == "NaN") sd.s <- 2 * s / 3
  
  if (skip.number == "NaN") skip.number <- 5
  
  number.iterations <- (skip.number + 1) * chain.length + burnin
  
  # run the MCMC chain.
  chain.output <- MCMC.chain(
    X, beta, g, s, i, p1, p2, lambda1, lambda2, sd.g, sd.s, skip.number,
    number.iterations, burnin, J.numb
  )
 
  return(chain.output)
  
  chain <- chain.output[, 1:3]
  acceptance.rate <- mean(chain.output[, 4])
  cat("acceptance rate: ", acceptance.rate)
  stop()
  # estimated threshold
  median.index <- median(chain.output[, 3])
  
  beta.length <- length(beta)
  
  beta.label <- paste("Beta", 1:beta.length)
  names(beta) <- beta.label
  
  # sort the quantiles chain.
  quantile.chain <- apply(cbind(chain.output[, 5:(4 + beta.length)], 0), 2, sort)
  
  q.est.median <- apply(quantile.chain, 2, median)[1:beta.length]
  names(q.est.median) <- beta.label
  
  quantile.chain <- as.matrix(quantile.chain[, 1:beta.length])
  
  gamma.est.median <- median(chain[, 1])
  
  
  # CI for the quantiles. Upper, lower, and symmetric
  quantile.CI.upper <- 
    quantile.chain[(1 - (1 - CI.level) / 2) * (chain.length), 1:beta.length]
  quantile.CI.lower <- 
    quantile.chain[(1 - CI.level) / 2 * (chain.length), 1:beta.length]
  quantile.CI.symmetric <- rbind(quantile.CI.lower, quantile.CI.upper)
  
  colnames(quantile.CI.symmetric) <- beta.label
  
  
  # CI for the "shortest" interval.
  quantile.CI.shortest <- matrix(0, 2, beta.length)
  for (i in 1:beta.length) {
    quantile.CI.shortest[, i] <- CI.short.fast(quantile.chain[, i], CI.level)
  }
  
  rownames(quantile.CI.shortest) <- 
    c("quintile.CI.short.lower", "quintile.CI.short.upper")
  colnames(quantile.CI.shortest) <- beta.label
  
  return(list(
    median.threshold = X[median.index], Beta = beta, Beta.quantile = q.est.median,
    quantile.CI.symmetric = quantile.CI.symmetric,
    quantile.CI.shortest = quantile.CI.shortest,
    gamma = gamma.est.median
  ))
}






#################################################################################################
#### These are the necessary functions to run the Fiducial.GPD.unknown.thresh function above######
#################################################################################################


###### Calcualte the Jacobian ##########

Jacobian <- function(g, s, a, J.numb, X) {
  # data must be trucated for X>a
  X.choose.3 <- matrix(0, J.numb, 3)
  for (i in 1:J.numb) {
    X.choose.3[i, ] <- sample(X, 3)
  }
  
  if (g != 0) {
    X.diff <- cbind(X.choose.3[, 2] - X.choose.3[, 3], -X.choose.3[, 1] +
                      X.choose.3[, 3], X.choose.3[, 1] - X.choose.3[, 2])
    J.mat <- (log(1 + g * (X.choose.3 - a) / s) * 
                (1 + g * (X.choose.3 - a) / s) / g^2) * (X.diff)
    J.mean <- mean(abs(apply(J.mat, 1, sum)))
  } else {
    X.diff <- cbind(X.choose.3[, 2] - X.choose.3[, 3], -X.choose.3[, 1] +
                      X.choose.3[, 3], X.choose.3[, 1] - X.choose.3[, 2])
    J.mat <- apply(X.diff, 1, prod) / (2 * s)
    J.mean <- mean(abs(J.mat))
  }
  
  return(J.mean)
}

############# log of the fiducial density for (gamma,sigma,a)###############

log.gpd.dens <- function(g, s, a, i, X, J.numb, n) {
  # data must be pre-sorted
  
  X <- X[X > a]
  
  if (s > 0 & g > (-s / max(X - a)) & min(X - a) > 0 & a > 0 & g > -.5) {
    J <- Jacobian(g, s, a, J.numb, X)
    if (g != 0) {
      log.density <- 
        sum((-1 / g - 1) * log(1 + g * (X - a) / s)) + log(J) - n * log(s + a)
    } else {
      log.density <- -1 / s * sum(X - a) + log(J) - n * log(s + a)
    }
  } else {
    log.density <- -Inf
  }
  return(log.density)
}


####### calcualte the Beta-quantile for a value or a vector of Beta values###########

quantile.value <- function(g, s, a, prob.a, beta) {
  alpha <- 1 - beta
  alpha <- (alpha / prob.a)
  if (g != 0) {
    Q <- a + s / g * ((alpha)^(-g) - 1)
  }
  if (g == 0) {
    Q <- a - s * log(alpha)
  }
  return(Q)
}

######## propose a new (gamma,sigma) value or a new index for the thershold##########

MCMC.newpoint <- function(g, s, i, p1, p2, lambda, sd.g, sd.s, X, J.numb, n) {
  # data must be sorted
  # p1 is the probability that we change g ans s.
  # p2 is the probability that we move + or - when we move indexes
  # lambda is the poisson parameter
  
  # this portion proposes a new i (new threshold)
  a <- X[i]
  if (runif(1) > p1) {
    Bern <- sample(c(-1, 1), 1, prob = c(1 - p2, p2))
    
    if (Bern == 1) {
      plus.minus <- n
      while (plus.minus > (n - i - 10)) {
        plus.minus <- rpois(1, lambda)
      }
      i.star <- i + plus.minus
      
      dens.pois.star <- 1 / ppois(n - i - 10, lambda) * p2
      dens.pois <- 1 / ppois(i.star - 1, lambda) * (1 - p2)
    }
    
    if (Bern == -1) {
      plus.minus <- n
      lambda <- min(i, lambda)
      while (plus.minus > (i - 1)) {
        plus.minus <- rpois(1, lambda)
      }
      i.star <- i - plus.minus
      dens.pois.star <- 1 / ppois(i - 1, lambda) * (1 - p2)
      dens.pois <- 1 / ppois(n - i.star - 10, lambda) * p2
    }
    
    a.star <- X[i.star]
    
    gs.star <- c(g, s + g * (X[i.star] - X[i]))
    MH.ratio <- 
      exp(log.gpd.dens(gs.star[1], gs.star[2], a.star, i.star, X, J.numb, n) - 
            log.gpd.dens(g, s, a, i, X, J.numb, n)) *
      dens.pois / dens.pois.star
    
    
    # this proposes a new gamma and sigma.
  } else {
    i.star <- i
    a.star <- X[i.star]
    gs.star <- rcauchy(2, c(g, s), c(sd.g, sd.s)) # rnorm(2,c(g,s),c(sd.g,sd.s))
    
    MH.ratio <- 
      exp(log.gpd.dens(gs.star[1], gs.star[2], a.star, i.star, X, J.numb, n) - 
            log.gpd.dens(g, s, a, i, X, J.numb, n))
  }
  
  U <- runif(1)
  
  # Metropolis-Hastings ratio to accept the new point.
  if (U <= MH.ratio & MH.ratio != "NaN" & MH.ratio != "Inf") {
    newpoint <- c(gs.star, i.star, 0)
  } else {
    newpoint <- c(g, s, i, 0)
  }
  
  return(newpoint)
}


########### The function that runs the MCMC chain we will use to create intervals and estimates##############

MCMC.chain <- function(X, beta, g, s, i, p1, p2, lambda1, lambda2, sd.g, sd.s, 
                       skip.number, number.iterations, burnin, J.numb) {
  n <- length(X)
  X <- sort(X)
  # transfrom the data X'=X-X(1)
  min.X <- min(X)
  X <- (X - min.X)
  
  x.t <- matrix(0, nrow = number.iterations, ncol = 4 + length(beta))
  x.t[1, ] <- c(g, s, i, 0, quantile.value(g, s, X[i], 1 - i / n, beta))
  
  # propose new points and run the MCMC chain for a specified chain length
  for (j in 1:(number.iterations - 1)) {
    if (j %% 10 == 1) {
      lambda <- lambda2
    } else {
      lambda <- lambda1
    }
    x.t[j + 1, 1:4] <- MCMC.newpoint(
      x.t[j, 1], x.t[j, 2], x.t[j, 3], p1, p2, lambda, sd.g, sd.s, X, J.numb, n
    )
    x.t[j + 1, 5:(4 + length(beta))] <- quantile.value(
      x.t[j + 1, 1], x.t[j + 1, 2], X[x.t[j + 1, 3]], 1 - x.t[j + 1, 3] / n, beta
    )
  }
  
  # Delet the first values as the "burnin"
  x.t <- x.t[(burnin + 1):number.iterations, ]
  
  # Thin the chain by keeping every skip.number+1 iteration of the MCMC chain
  number.iterations <- nrow(x.t)
  every.ith <- c(1, rep(0, skip.number))
  eliminate.vector <- 
    rep(every.ith, 
        ceiling((number.iterations) / length(every.ith)))[1:number.iterations]
  
  x.t <- x.t[eliminate.vector == 1, ]
  
  # Indictor for the acceptance rate.
  for (i in 1:(length(x.t[, 1]) - 1)) {
    if (x.t[i, 5] != x.t[i + 1, 5]) {
      x.t[i + 1, 4] <- 1
    }
  }
  
  x.t[, 5:(4 + length(beta))] <- x.t[, 5:(4 + length(beta))] + min.X
  
  return(x.t)
}


###### CI for the "shortest" interval.############

CI.short.fast <- function(chain, confidence.level) {
  chain.length <- length(chain)
  chain <- sort(chain)
  
  num.short.int <- min(chain.length * (1 - confidence.level) + 1, chain.length)
  length.short.int <- rep(0, num.short.int)
  
  length.short.int <- chain[(chain.length - num.short.int + 1):chain.length] - chain[1:num.short.int]
  index <- rank(length.short.int)
  index <- (index == min(index)) * 1:num.short.int
  lower.ix <- floor(median(index[index != 0]))
  shortest.int <- c(lower = chain[lower.ix], upper = chain[chain.length - num.short.int + lower.ix])
  
  return(shortest.int)
}


data(rain)
set.seed(666)
X <- rgamma(10000, shape = 4, scale = 1)  #rain#rgamma(1000, 1, 1)

xt <- Fiducial.GPD.unknown.thresh(
  X, c(0.25,0.5,0.75), 0.95, chain.length = 10000, burnin = 2000
)

coda <- as.mcmc(xt)
summary(coda)
qgamma(c(0.25,0.5,0.75),4)
