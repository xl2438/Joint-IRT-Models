######################### Setting item parameters #####################
n_item <- 15
a <- rlnorm(n_item, 0, 0.15)
b <- runif(n_item, -2, 2.0)
ga <- rnorm(n_item, 0, 1)
alpha <- rlnorm(n_item, 0, 0.15)
beta <- rnorm(n_item, 0, 1)
#################################################################
one_rep <- function(rep, n_person, a, b, ga, alpha, beta, corr) {
  p_irt <- function(theta, a, b, ga, w) {
    out <- 1 / (1 + exp(a * (outer(b, theta, "-") - ga * w)))
    return(out)
  }
  irt_2pl <- function(theta, a, b) {
    out <- 1 / (1 + exp(a * outer(b, theta, "-")))
    return(out)
  }
  Rcpp::sourceCpp("test.cpp")
  n_item <- length(a)
  rho <- matrix(c(1, corr, corr, 1), nrow = 2)
  rho_sd <- c(1, 1)
  sig <- sweep(sweep(rho, 1L, rho_sd, "*"), 2L, rho_sd, "*")
  theta <- mvtnorm::rmvnorm(n_person, mean = c(0, 0), sigma = sig)
  w <- 1 * (irt_2pl(theta[, 2], alpha, beta) > runif(n_item * n_person, 0, 1))
  x <- 1 * (p_irt(theta[, 1], a, b, ga, w) > runif(n_item * n_person, 0, 1))
  s <- estimate(x, w)
  return(s)
}
################################################
n_rep <- 10
corr <- 0.3
for (i in 1:30) {
  temp <- one_rep(1, 1000, a, b, ga, alpha, beta, corr)
  if (any(abs(temp$a) > 2)) {
    break
  }
  print(i)
}

x <- temp$x
w <- temp$w
s <- estimate(x, w)

sum(x[6, ])
sum(w[6, ])

table(x[6, ], w[6, ])