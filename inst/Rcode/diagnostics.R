library(gfiExtremes)
library(ggmcmc)

set.seed(1111111111L)
X <- rweibull(100L, shape = 1.5)
summary(X)
plot(X, pch = 19)

gf <- gfigpd2(
  X, beta = c(99, 99.5, 99.9)/100,# threshold.init = 0.7,
  allParameters = TRUE,
  burnin = 20000L, iter = 20000L, thin = 10L
)
summary(gf)
qweibull(c(99, 99.5, 99.9)/100, shape = 1.5)

pChains <- gf[, c("gamma", "sigma")]
qChains <- as.mcmc.list(
  lapply(gf[, c("beta1", "beta2", "beta3")], log)
)

pgg <- ggs(pChains)
qgg <- ggs(qChains)

ggs_density(pgg)
ggs_density(qgg)

ggs_traceplot(pgg)
ggs_traceplot(qgg)

ggs_running(pgg)
ggs_running(qgg)

ggs_compare_partial(pgg)
ggs_compare_partial(qgg)

ggs_autocorrelation(pgg)
ggs_autocorrelation(qgg)

ggs_Rhat(pgg) + xlab("R_hat")
ggs_Rhat(qgg) + xlab("R_hat")

ggs_geweke(qgg)

heidel.diag(qChains)

