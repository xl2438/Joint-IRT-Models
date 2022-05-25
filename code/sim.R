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
library(parallel)
n_rep <- 100
start_time <- Sys.time()
cl <- makeCluster(10)
corr <- 0.0
out_800 <- clusterApply(cl, c(1:n_rep), one_rep, n_person = 800, a = a, b = b,
 ga = ga, alpha = alpha, beta = beta, corr = corr)
out_1500 <- clusterApply(cl, c(1:n_rep), one_rep, n_person = 1500, a = a, b = b,
 ga = ga, alpha = alpha, beta = beta, corr = corr)
out_3000 <- clusterApply(cl, c(1:n_rep), one_rep, n_person = 3000, a = a, b = b,
 ga = ga, alpha = alpha, beta = beta, corr = corr)
save.image(file = "cor0.RData")
corr <- 0.3
out_800 <- clusterApply(cl, c(1:n_rep), one_rep, n_person = 800, a = a, b = b,
 ga = ga, alpha = alpha, beta = beta, corr = corr)
out_1500 <- clusterApply(cl, c(1:n_rep), one_rep, n_person = 1500, a = a, b = b,
 ga = ga, alpha = alpha, beta = beta, corr = corr)
out_3000 <- clusterApply(cl, c(1:n_rep), one_rep, n_person = 3000, a = a, b = b,
 ga = ga, alpha = alpha, beta = beta, corr = corr)
save.image(file = "cor3.RData")
corr <- 0.6
out_800 <- clusterApply(cl, c(1:n_rep), one_rep, n_person = 800, a = a, b = b,
 ga = ga, alpha = alpha, beta = beta, corr = corr)
out_1500 <- clusterApply(cl, c(1:n_rep), one_rep, n_person = 1500, a = a, b = b,
 ga = ga, alpha = alpha, beta = beta, corr = corr)
out_3000 <- clusterApply(cl, c(1:n_rep), one_rep, n_person = 3000, a = a, b = b,
 ga = ga, alpha = alpha, beta = beta, corr = corr)
save.image(file = "cor6.RData")
stopCluster(cl)
end_time <- Sys.time()
end_time - start_time

##################################################
a_out <- array(NA, dim = c(15, 3, 3, n_rep)) #item, cor condition, size, rep
b_out <- array(NA, dim = c(15, 3, 3, n_rep))
alpha_out <- array(NA, dim = c(15, 3, 3, n_rep))
beta_out <- array(NA, dim = c(15, 3, 3, n_rep))
ga_out <- array(NA, dim = c(15, 3, 3, n_rep))
rho_out <- array(NA, dim = c(1, 3, 3, n_rep))
load("cor0.RData")
for (i in 1:100) {
  a_out[, 1, 1, i] <- out_800[[i]]$a
  a_out[, 1, 2, i] <- out_1500[[i]]$a
  a_out[, 1, 3, i] <- out_3000[[i]]$a
  b_out[, 1, 1, i] <- out_800[[i]]$b
  b_out[, 1, 2, i] <- out_1500[[i]]$b
  b_out[, 1, 3, i] <- out_3000[[i]]$b
  alpha_out[, 1, 1, i] <- out_800[[i]]$alpha
  alpha_out[, 1, 2, i] <- out_1500[[i]]$alpha
  alpha_out[, 1, 3, i] <- out_3000[[i]]$alpha
  beta_out[, 1, 1, i] <- out_800[[i]]$beta
  beta_out[, 1, 2, i] <- out_1500[[i]]$beta
  beta_out[, 1, 3, i] <- out_3000[[i]]$beta
  ga_out[, 1, 1, i] <- out_800[[i]]$ga
  ga_out[, 1, 2, i] <- out_1500[[i]]$ga
  ga_out[, 1, 3, i] <- out_3000[[i]]$ga
  rho_out[, 1, 1, i] <- out_800[[i]]$rho
  rho_out[, 1, 2, i] <- out_1500[[i]]$rho
  rho_out[, 1, 3, i] <- out_3000[[i]]$rho
}
load("cor3.RData")
for (i in 1:100) {
  a_out[, 2, 1, i] <- out_800[[i]]$a
  a_out[, 2, 2, i] <- out_1500[[i]]$a
  a_out[, 2, 3, i] <- out_3000[[i]]$a
  b_out[, 2, 1, i] <- out_800[[i]]$b
  b_out[, 2, 2, i] <- out_1500[[i]]$b
  b_out[, 2, 3, i] <- out_3000[[i]]$b
  alpha_out[, 2, 1, i] <- out_800[[i]]$alpha
  alpha_out[, 2, 2, i] <- out_1500[[i]]$alpha
  alpha_out[, 2, 3, i] <- out_3000[[i]]$alpha
  beta_out[, 2, 1, i] <- out_800[[i]]$beta
  beta_out[, 2, 2, i] <- out_1500[[i]]$beta
  beta_out[, 2, 3, i] <- out_3000[[i]]$beta
  ga_out[, 2, 1, i] <- out_800[[i]]$ga
  ga_out[, 2, 2, i] <- out_1500[[i]]$ga
  ga_out[, 2, 3, i] <- out_3000[[i]]$ga
  rho_out[, 2, 1, i] <- out_800[[i]]$rho
  rho_out[, 2, 2, i] <- out_1500[[i]]$rho
  rho_out[, 2, 3, i] <- out_3000[[i]]$rho
}
load("cor6.RData")
for (i in 1:100) {
  a_out[, 3, 1, i] <- out_800[[i]]$a
  a_out[, 3, 2, i] <- out_1500[[i]]$a
  a_out[, 3, 3, i] <- out_3000[[i]]$a
  b_out[, 3, 1, i] <- out_800[[i]]$b
  b_out[, 3, 2, i] <- out_1500[[i]]$b
  b_out[, 3, 3, i] <- out_3000[[i]]$b
  alpha_out[, 3, 1, i] <- out_800[[i]]$alpha
  alpha_out[, 3, 2, i] <- out_1500[[i]]$alpha
  alpha_out[, 3, 3, i] <- out_3000[[i]]$alpha
  beta_out[, 3, 1, i] <- out_800[[i]]$beta
  beta_out[, 3, 2, i] <- out_1500[[i]]$beta
  beta_out[, 3, 3, i] <- out_3000[[i]]$beta
  ga_out[, 3, 1, i] <- out_800[[i]]$ga
  ga_out[, 3, 2, i] <- out_1500[[i]]$ga
  ga_out[, 3, 3, i] <- out_3000[[i]]$ga
  rho_out[, 3, 1, i] <- out_800[[i]]$rho
  rho_out[, 3, 2, i] <- out_1500[[i]]$rho
  rho_out[, 3, 3, i] <- out_3000[[i]]$rho
}

a_mean <- apply(a_out, c(1, 2, 3), mean)
b_mean <- apply(b_out, c(1, 2, 3), mean)

out <- matrix(unlist(out_1500), nrow = n_item * 5 + 1)
a_out <- out[1:15, ]
apply(a_out, 1, mean)