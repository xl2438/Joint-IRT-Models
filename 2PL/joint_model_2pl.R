n_item <- 20
n_person <- 3000

rho <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
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

compute_joint_dist <- function(ind, q_pts, q_w, x, w, a, b, alpha, beta, ga) {
  p_x <- p_irt(rep(q_pts[ind, 1], dim(w)[2]), a, b, ga, w)
  p_w <- irt_2pl(rep(q_pts[ind, 2], dim(w)[2]), alpha, beta)
  l_x <- x * log(p_x) + (1 - x) * log(1 - p_x)
  l_w <- w * log(p_w) + (1 - w) * log(1 - p_w)
  l <- apply(l_x, 2, sum) + apply(l_w, 2, sum) + q_w[ind]
  return(exp(l))
}
compute_p_w <- function(x, w, quad_pts, quad_w, a, b, alpha, beta, ga) {
  post_w <- sapply(seq_len(length(quad_w)), FUN = compute_joint_dist, q_pts =
   quad_pts, q_w = quad_w, x = x, w = w, a = a, b = b, alpha = alpha, beta =
   beta, ga = ga)
  z <- apply(post_w, 1, sum)
  post_w <- post_w / z
  return(t(post_w))
}
compute_res <- function(q_pts, a, b, x) {
  p <- irt_2pl(q_pts, a, b)
  out <- outer(as.vector(x), as.vector(p), "-")
  return(t(out))
}
obj_fun_j <- function(v, x_vec, q_pts, post_w) {
  a <- v[1]
  b <- v[2]
  out <- rep(NA, 2)
  res <- compute_res(q_pts, a, b, x_vec)
  out[1] <- sum(res * post_w)
  out[2] <- sum((q_pts - b) * res * post_w)
  return(out)
}
jacobian_obj <- function(v, x_vec, q_pts, post_w) {
  a <- v[1]
  b <- v[2]
  x_vec <- as.vector(x_vec)
  out <- matrix(NA, nrow = 2, ncol = 2)
  p <- as.vector(irt_2pl(q_pts, a, b))
  pq <- p * (1 - p)
  res <- q_pts - b
  out[1, 1] <- -sum(res * pq * post_w)
  out[1, 2] <- sum(a * pq * post_w)
  out[2, 1] <- -sum(res^2 * pq * post_w)
  out[2, 2] <- sum((outer(p, x_vec, "-") +  a * pq * res) * post_w)
  return(out)
}
estimate_par_j <- function(ind, x, q_pts, post_w) {
  sln <- nleqslv::nleqslv(x = c(1, 0), fn = obj_fun_j, jac = jacobian_obj, x_vec
   = x[ind, ], q_pts = q_pts, post_w = post_w, method = "Newton")$x
  return(sln)
}
estimate_item_pars <- function(x, q_pts, post_w) {
  out <- sapply(c(1:dim(x)[1]), FUN = estimate_par_j, x = x, q_pts = q_pts,
    post_w = post_w)
  return(out)
}
#############################################################
compute_res_ga <- function(q_pts, a, b, ga, x, w) {
  p <- 1 / (1 + exp(-a * outer(q_pts - b, w * ga, "+")))
  out <- x - t(p)
  return(t(out))
}
obj_fun_ga_j <- function(v, x_vec, w_vec, q_pts, post_w) {
  a <- v[1]
  b <- v[2]
  ga <- v[3]
  out <- rep(NA, 3)
  res <- compute_res_ga(q_pts, a, b, ga, x_vec, w_vec)
  out[1] <- sum(res * post_w)
  out[2] <- sum(outer(q_pts - b, w_vec * ga, "+") * res * post_w)
  out[3] <- sum(t(w_vec * t(res)) * post_w)
  return(out)
}
jacobian_obj_ga_j <- function(v, x_vec, w_vec, q_pts, post_w) {
  a <- v[1]
  b <- v[2]
  ga <- v[3]
  out <- matrix(NA, nrow = 3, ncol = 3)
  temp_logit <- outer(q_pts - b, w_vec * ga, "+")
  p <- sapply(w_vec, FUN = p_irt, theta = q_pts, a = a, b = b, ga = ga)
  p_q <- p * (1 - p)
  out[1, 1] <- -sum(temp_logit * p_q * post_w)
  out[1, 2] <- sum(a * p_q * post_w)
  out[1, 3] <- -sum(a * w_vec * t(p_q * post_w))
  out[2, 1] <- -sum(temp_logit^2 * p_q * post_w)
  out[2, 2] <- sum((t(t(p) - x_vec) + temp_logit * a * p_q) * post_w)
  out[2, 3] <- sum(t(w_vec * (x_vec - t(p)) - w_vec * t(temp_logit * a * p_q)) *
   post_w)
  out[3, 1] <- -sum(w_vec * t(temp_logit * p_q * post_w))
  out[3, 2] <- -out[1, 3]
  out[3, 3] <- out[1, 3]
  return(out)
}
estimate_par_ga_j <- function(ind, x, w, q_pts, post_w) {
  sln <- nleqslv::nleqslv(x = c(1, 0, 0), fn = obj_fun_ga_j, jac =
   jacobian_obj_ga_j, x_vec = as.vector(x[ind, ]), w_vec = as.vector(w[ind, ]),
   q_pts = q_pts, post_w = post_w, method = "Newton")$x
  return(sln)
}
estimate_item_pars_ga <- function(x, w, q_pts, post_w) {
  out <- sapply(c(1:dim(x)[1]), FUN = estimate_par_ga_j, x = x, w = w, q_pts =
    q_pts, post_w = post_w)
  return(out)
}
obj_rho <- function(rho, expec_theta_squared, expec_theta_prod, n) {
  out <- n * rho^3 - expec_theta_prod * rho^2 - (1 - expec_theta_squared / n) *
  n * rho - expec_theta_prod
  return(out)
}
estimate_rho <- function(quad_pts, post_w) {
  n_person <- dim(post_w)[2]
  w_q <- apply(post_w, 1, sum)
  expec_theta_squared <- sum(apply(quad_pts^2, 1, sum) * w_q)
  expec_theta_prod <- sum(apply(quad_pts, 1, prod) * w_q)
  rho <- uniroot(f = obj_rho, interval = c(-1, 1), expec_theta_squared =
   expec_theta_squared, expec_theta_prod = expec_theta_prod, n = n_person)$root
  return(rho)
}
################## EM initialization ########################################
a_est <- rlnorm(n_item, 0, 0.15)
b_est <- rnorm(n_item)
ga_est <- rnorm(n_item)
alpha_est <- rlnorm(n_item, 0, 0.15)
beta_est <- rnorm(n_item)
rho_est <- 0.4
rho_sd <- c(1, 1)
##############################################################################

for (i in 1:15) {
  sig <- sweep(sweep(matrix(c(1, rho_est, rho_est, 1), nrow = 2), 1L, rho_sd,
   "*"), 2L, rho_sd, "*")
  temp <- MultiGHQuad::init.quad(Q = 2, prior = list(mu = c(0, 0), Sigma =
    matrix(c(1, rho_est, rho_est, 1), nrow = 2, ncol = 2)), ip = 13)
  quad_w <- temp$W #log weights
  quad_pts <- temp$X
  post_w <- compute_p_w(x = x, w = w, quad_pts = quad_pts, quad_w = quad_w, a =
    a, b = b, alpha = alpha, beta = beta, ga = ga)
  # post_w <- compute_posterior(quad_pts, quad_w, x, w, a_est, b_est, alpha_est,
  #  beta_est, ga_est)
  temp <- estimate_item_pars(x = w, q_pts = as.vector(quad_pts[, 2]), post_w =
   post_w)
  alpha_est <- temp[1, ]
  beta_est <- temp[2, ]
  temp <- estimate_item_pars_ga(x, w, quad_pts[, 1], post_w)
  a_est <- temp[1, ]
  b_est <- temp[2, ]
  ga_est <- temp[3, ]
  rho_est <- estimate_rho(quad_pts, post_w)
  print(i)
}

Rcpp::sourceCpp("compute_posterior_weights.cpp")
temp <- compute_quad_pts(0.4, 13)

temp <- update_quad_weights(0.4)
str(temp)



temp <- compute_posterior(quad_pts, quad_w, x, w, a, b, alpha, beta, ga)

j <- 12
estimate_a_b_gamma(x[j, ], w[j, ], quad_pts[, 1], temp)


update_a_b_gamma(x, w, quad_pts, temp)
cbind(alpha, beta)
cbind(a, b, ga)

estimate_a_b_gamma(x[20, ], w[20, ], quad_pts[, 1], temp)
estimate_par_ga_j(8, x, w, quad_pts[, 1], temp)


start_time <- Sys.time()
update_a_b_gamma(x, w, quad_pts[, 1], temp)
Sys.time() - start_time

start_time <- Sys.time()
for(i in 1:1000) tmp1 <- obj_fun_ga_j(c(a[1], b[1], ga[1]), x[1, ], w[1, ], quad_pts[,1], temp)
Sys.time() - start_time

tmp1 <- obj_fun_ga_j(c(a[1], b[1], ga[1]), x[1, ], w[1, ], quad_pts[,1], temp)



# start_time <- Sys.time()
# tmp <- update_alpha_beta(w, quad_pts[, 2], temp)
# Sys.time() - start_time

# sum(abs(tmp - cbind(alpha, beta)))

# start_time <- Sys.time()
# for (i in 1:200000){
#   tmp <- irt_2pl(theta[,2], alpha, beta)  
# }
# Sys.time() - start_time
