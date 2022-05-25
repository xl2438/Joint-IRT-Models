n_item <- 20
n_person <- 2000

rho <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
rho_sd <- c(1, 1)
sig <- sweep(sweep(rho, 1L, rho_sd, "*"), 2L, rho_sd, "*")

theta <- mvtnorm::rmvnorm(n_person, mean = c(0, 0), sigma = sig)
b <- runif(n_item, -2, 2.0)
ga <- rnorm(n_item, 0, 1)
beta <- rnorm(n_item, 0, 1)

irt <- function(theta, b, beta, w) {
  out <- matrix(NA, nrow = length(theta), ncol = length(b))
  for (i in 1:dim(out)[1]) {
    for (j in 1:dim(out)[2]) {
      out[i, j] <- theta[i] - b[j] + w[i, j] * beta[j]
    }
  }
  p <- 1 / (1 + exp(-1 * out))
  return(p)
}

cal <- function(tau, ga) {
  out <- 1 / (1 + exp(outer(ga, tau, "-")))
  return(t(out))
}

w <- 1 * (cal(theta[, 2], ga) > runif(n_item * n_person, 0, 1))
x <- 1 * (irt(theta[, 1], b, beta, w) > runif(n_item * n_person, 0, 1))

dim(w)
dim(x)
##########################################
dat <- list(N = n_person, J = n_item, x = x, w = w)
library(rstan)
options(mc.cores = parallel::detectCores())
fit <- stan("tool_usage.stan", data = dat, chains = 2, iter = 4000)
print(fit, pars = c("rho", "b", "beta", "ga"))
temp <- extract(fit, pars = c("rho", "b", "beta", "ga"));
b_eap <- apply(temp[[2]], 2, mean)
b_eap
beta_eap <- apply(temp[[3]], 2, mean)
beta_eap
beta
