library(ismev)
library(coda)

######## This is the actual function to run for the package############


test <- function(
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
  
  X <- sort(X) # -->> so there's no need to sort in C++
  n <- length(X)
  
  # Intialize the default values for the tuning parameters of the MCMC chain.
  
  if (i == "NaN") i <- floor(n * .85)
  
  if(is.nan(g) || is.nan(s)){
    mle.fit <- gpd.fit(X, X[i], show = FALSE)
    if (is.nan(g)) g <- mle.fit$mle[2L]
    if (is.nan(s)) s <- mle.fit$mle[1L]
  }
  
  if (sd.g == "NaN") sd.g <- 2 * abs(g) / 3
  if (sd.s == "NaN") sd.s <- 2 * s / 3
  
  if (skip.number == "NaN") skip.number <- 5L
  
  number.iterations <- (skip.number + 1L) * chain.length + burnin
  
  # run the MCMC chain.
  x.t <- gfconfExtremes:::MCMCchain(
    X, beta, g, s, i, p1, p2, lambda1, lambda2, sd.g, sd.s, skip.number,
    number.iterations, burnin, J.numb, 666L 
  ) # skip.number and burnin not used
  
  # Thin the chain by keeping every skip.number+1 iteration of the MCMC chain
  number.iterations <- nrow(x.t)
  every.ith <- c(TRUE, rep(FALSE, skip.number))
  eliminate.vector <- 
    rep(every.ith, 
        ceiling((number.iterations) / length(every.ith)))[1:number.iterations]
  
  x.t <- x.t[eliminate.vector, ]
  
  print(dim(x.t))
  
  # Indictor for the acceptance rate.
  ii <- integer(nrow(x.t))
  for (i in 1L:(nrow(x.t) - 1L)) {
    if (x.t[i, 4L] != x.t[i + 1L, 4L]) {
      ii[i + 1L] <- 1L
    }
  }
  
  x.t[, 4L:(3L + length(beta))] <- x.t[, 4L:(3L + length(beta))] + X[1L]
  
  acceptance.rate <- mean(ii)
  cat("acceptance rate: ", acceptance.rate)
  
  
  return(x.t)
  
}

data("rain")
set.seed(666)
U <- runif(50000,0,10)
P <- 10 + 1*(runif(50000)^(-0.5)-1)/0.5#Pareto::rPareto(10000, 10, 1)  #rain#rgamma(1000, 1, 1)

set.seed(666)
a <- 10
U <- runif(50000)
xi <- 0.5
P <- a + 1*(U^(-xi)-1)/xi#Pareto::rPareto(10000, 10, 1)  #rain#rgamma(1000, 1, 1)

# p * 1/a = 1 - p => p = 1/(1/a+1)
p <- 1/(1/a+1)
X <- numeric(50000)
for(i in 1:50000){
  X[i] <- sample(c(a*U[i],P[i]),1, prob = c(p,1-p))
}
system.time(xt <- test(X, c(0.9, 0.95, 0.99), chain.length = 10000, burnin = 2000))
#str(X)

coda <- as.mcmc(xt)
summary(coda)
#Pareto::qPareto(c(0.25,0.5,0.75), 10, 1)
quantile(X, c(0.9, 0.95, 0.99))