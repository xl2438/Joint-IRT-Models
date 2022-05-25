n_item <- 20
n_person <- 3000

rho <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
rho_sd <- c(1, 1)
sig <- sweep(sweep(rho, 1L, rho_sd, "*"), 2L, rho_sd, "*")

theta <- mvtnorm::rmvnorm(n_person, mean = c(0, 0), sigma = sig)
a <- rlnorm(n_item, 0, 0.15)
b <- runif(n_item, -2, 2.0)
ga <- rnorm(n_item, 0, 1)
alpha <- rlnorm(n_item, 0, 0.15)
beta <- rnorm(n_item, 0, 1)
###### compute posterior of latent variabels ######
p_irt <- function(theta, a, b, ga, w) {
  out <- 1 / (1 + exp(a * (outer(b, theta, "-") - ga * w)))
  return(out)
}
irt_2pl <- function(theta, a, b) {
  out <- 1 / (1 + exp(a * outer(b, theta, "-")))
  return(out)
}
w <- 1 * (irt_2pl(theta[, 2], alpha, beta) > runif(n_item * n_person, 0, 1))
x <- 1 * (p_irt(theta[, 1], a, b, ga, w) > runif(n_item * n_person, 0, 1))

library(Rcpp)
sourceCpp("test.cpp")

fit <- estimate(x, w) #default to max iterations = 40 or the relative
# improvement in loglikelihood is less than 1e-6.

fit$rho
fit$beta