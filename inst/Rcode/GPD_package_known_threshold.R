library(ismev)

######## This is the actual function to run for the package############

Fiducial.GPD.known.thresh <- function(
                                      X, a, beta, CI.level, g = NaN, s = NaN, sd.g = NaN,
                                      sd.s = NaN, skip.number = NaN, chain.length = 10000,
                                      burnin = 2000, J.numb = 50) {
  # X is the vector of data.
  # beta is the value or vector of to estimate the Beta-quantile.
  # CI.level is the confidence level for confidence intervals.
  # sd.g is the standard deviation for the proposed gamma.
  # sd.s is the standard deviation for the proposed sigma.
  # skip.number is a thinning number for the MCMC chain. (e.g. if skip.number is 1 we will use ever other MCMC iteration).
  # chain.length is the length of the MCMC chain.
  # burnin is number of the first MCMC iterations omitted from the chain.
  # J.numb is the number of subsamples that are taken from the Jacobian.

  X <- sort(X)



  # Intialize the default values for the tuning parameters of the MCMC chain.

  mle.fit <- gpd.fit(X, a, show = FALSE)

  n <- sum(X >= a)
  if (g == "NaN") g <- mle.fit$mle[2]
  if (s == "NaN") s <- mle.fit$mle[1]

  if (sd.g == "NaN") sd.g <- .3 / log(n, 20)
  if (sd.s == "NaN") sd.s <- s * .3 / log(n, 20)

  if (skip.number == "NaN") skip.number <- 1


  num.iterations <- (skip.number + 1) * chain.length + burnin


  # run the MCMC chain.
  prob.greater.a <- mean(X > a)

  chain.output <- MCMC.chain(
    g, s, a, sd.g, sd.s, skip.number,
    num.iterations, beta, burnin, X, prob.greater.a, J.numb
  )

  chain <- chain.output[, 1:2]
  acceptance.rate <- mean(chain.output[, 3])

  beta.length <- length(beta)

  beta.label <- paste("Beta", 1:beta.length)
  names(beta) <- beta.label

  # sort the quantiles chain.
  quantile.chain <- apply(cbind(chain.output[, 4:(3 + beta.length)], 0), 2, sort)

  q.est.median <- apply(quantile.chain, 2, median)[1:beta.length]
  names(q.est.median) <- beta.label

  quantile.chain <- as.matrix(quantile.chain[, 1:beta.length])

  gamma.est.median <- median(chain[, 1])


  # CI for the quantiles. Upper, lower, and symmetric
  quantile.CI.upper <- quantile.chain[(1 - (1 - CI.level) / 2) * (chain.length), 1:beta.length]
  quantile.CI.lower <- quantile.chain[(1 - CI.level) / 2 * (chain.length), 1:beta.length]
  quantile.CI.symmetric <- rbind(quantile.CI.lower, quantile.CI.upper)

  colnames(quantile.CI.symmetric) <- beta.label


  # CI for the "shortest" interval.
  quantile.CI.shortest <- matrix(0, 2, beta.length)
  for (i in 1:beta.length) {
    quantile.CI.shortest[, i] <- CI.short.fast(quantile.chain[, i], CI.level)
  }

  rownames(quantile.CI.shortest) <- c("quintile.CI.short.lower", "quintile.CI.short.upper")
  colnames(quantile.CI.shortest) <- beta.label

  return(list(
    Beta = beta, Beta.quantile = q.est.median,
    quantile.CI.symmetric = quantile.CI.symmetric,
    quantile.CI.shortest = quantile.CI.shortest,
    gamma = gamma.est.median, acceptance.rate
  ))
}








#################################################################################################
#### These are the necessary functions to run the Fiducial.GPD.unknown.thresh function above######
#################################################################################################



Jacobian <- function(g, s, a, J.numb, X, n) {
  # data must be trucated for X>a
  if (n >= 250) {
    X.choose.2 <- matrix(0, J.numb, 2)
    for (i in 1:J.numb) {
      X.choose.2[i, ] <- sample(X, 2)
    }
  }


  if (g != 0) {
    if (n >= 250) {
      J.mat <- ((X.choose.2[, 1] - a) * (1 + g * (X.choose.2[, 2] - a) / s) * log(1 + g * (X.choose.2[, 2] - a) / s) -
        (X.choose.2[, 2] - a) * (1 + g * (X.choose.2[, 1] - a) / s) * log(1 + g * (X.choose.2[, 1] - a) / s)) / g^2
      J.mean <- mean(abs(J.mat))
    } else {
      X <- X - a
      X <- cbind(X)


      Xi.Xj <- X %*% t((1 + g * X / s) * log(1 + g * X / s))
      Xi.Xj.transp <- t(Xi.Xj)

      upper.tri <- matrix(0, n, n)
      upper.tri[upper.tri(upper.tri)] <- 1

      J.mean <- choose(n, 2)^(-1) *
        sum(abs((Xi.Xj - Xi.Xj.transp) * upper.tri)) / g^2
    }
  }
  if (g == 0) {
    if (n >= 250) {
      J.mat <- (X.choose.2[, 1] - a) * (X.choose.2[, 2] - a) * (X.choose.2[, 1] - X.choose.2[, 2]) / (2 * s^2)
      J.mean <- mean(abs(J.mat))
    } else {
      X <- X - a
      Xi.Xj <- X^2 %*% t(X)
      Xi.Xj.transp <- t(Xi.Xj)

      upper.tri <- matrix(0, n, n)
      upper.tri[upper.tri(upper.tri)] <- 1

      J.mean <- choose(n, 2)^(-1) *
        sum(abs((Xi.Xj - Xi.Xj.transp) * upper.tri) / (2 * s^2))
    }
  }
  return(J.mean)
}


############# log of the fiducial density for (gamma,sigma)###############


log.gpd.dens <- function(g, s, a, X, J.numb) {
  # data must be pre-sorted

  X <- X[X > a]
  n <- length(X)

  if (s > 0 & g > (-s / max(X - a))) {
    if (g != 0) {
      J <- Jacobian(g, s, a, J.numb, X, n)
      log.density <- sum(-log(s) + (-1 / g - 1) * log(1 + g * (X - a) / s)) + log(J)
    }
    if (g == 0) {
      J <- Jacobian(g, s, a, J.numb, X, n)
      log.density <- -n * log(s) - sum(X - a) / s + log(J)
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

MCMC.newpoint <- function(g, s, a, sd.g, sd.s, X, J.numb) {
  # data must be sorted


  gs.star <- rcauchy(2, c(g, s), c(sd.g, sd.s)) # rnorm(2,c(g,s),c(sd.g,sd.s))
  MH.ratio <- exp(log.gpd.dens(gs.star[1], gs.star[2], a, X, J.numb) - log.gpd.dens(g, s, a, X, J.numb))

  U <- runif(1)

  if (U <= MH.ratio & MH.ratio != "NaN" & MH.ratio != "Inf") {
    newpoint <- c(gs.star, 0)
  } else {
    newpoint <- c(g, s, 0)
  }

  return(newpoint)
}

########### The function that runs the MCMC chain we will use to create intervals and estimates##############


MCMC.chain <- function(g, s, a, sd.g, sd.s, skip.number,
                       number.iterations, beta, burnin, X, prob.greater.a, 
                       J.numb) {
  x.t <- matrix(0, nrow = number.iterations, ncol = 3 + length(beta))
  x.t[1, ] <- c(g, s, 0, quantile.value(g, s, a, prob.greater.a, beta))


  for (j in 1:(number.iterations - 1)) {
    x.t[j + 1, 1:3] <- 
      MCMC.newpoint(x.t[j, 1], x.t[j, 2], a, sd.g, sd.s, X, J.numb)
    x.t[j + 1, 4:(3 + length(beta))] <- 
      quantile.value(x.t[j + 1, 1], x.t[j + 1, 2], a, prob.greater.a, beta)
    #    if(x.t[j+1,4]!=x.t[j,4]){x.t[j+1,3]<-1}
  }

  x.t <- x.t[(burnin + 1):number.iterations, ]

  number.iterations <- nrow(x.t)
  every.ith <- c(1, rep(0, skip.number))
  eliminate.vector <- rep(every.ith, ceiling((number.iterations) / length(every.ith)))[1:number.iterations]

  x.t <- x.t[eliminate.vector == 1, ]

  for (i in 1:(length(x.t[, 1]) - 1)) {
    if (x.t[i, 4] != x.t[i + 1, 4]) {
      x.t[i + 1, 3] <- 1
    }
  }

  return(x.t)
}



MCMC.lines <- function(values.matrix, int.direction, int.area, match.area) {

  ## values.matrix=[g.values,s.values,acceptance]
  ## int.direction='g' or 's' depending on if you want dg or ds
  ## int.side='low' or 'up' depinding if you want the lower or upper area
  ## match.area=min(init.area.g,init.area.s, 1-init.area.g, 1-init.area.s)

  if (int.direction == "g") {
    if (int.area == "low") {
      g.values <- sort(values.matrix[, 1])
      value <- g.values[max((match.area) * length(g.values), 1)]
    }

    if (int.area == "up") {
      g.values <- sort(values.matrix[, 1])
      value <- g.values[(1 - match.area) * length(g.values)]
    }
  }

  if (int.direction == "s") {
    if (int.area == "low") {
      s.values <- sort(values.matrix[, 2])
      value <- s.values[max((match.area) * length(s.values), 1)]
    }

    if (int.area == "up") {
      s.values <- sort(values.matrix[, 2])
      value <- s.values[(1 - match.area) * length(s.values)]
    }
  }
  return(value)
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
