n_item <- 20
n_person <- 1000

rho <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
rho_sd <- c(1, 1)
sig <- sweep(sweep(rho, 1L, rho_sd, "*"), 2L, rho_sd, "*")

theta <- mvtnorm::rmvnorm(n_person, mean = c(0, 0), sigma = sig)
b <- runif(n_item, -2, 2.0)
ga <- rnorm(n_item, 0, 1)
beta <- rnorm(n_item, 0, 1)
###### compute posterior of latent variabels ######
p_irt <- function(theta, b, beta, w) {
  out <- 1 / (1 + exp(outer(b, theta, "-") - beta * w))
  return(out)
}
p_cal <- function(theta, ga) {
  out <- 1 / (1 + exp(outer(ga, theta, "-")))
  return(out)
}
p_irt_no_w <- function(b, theta) {
  out <- 1 / (1 + exp(outer(b, theta, "-")))
  return(out)
}
w <- 1 * (p_cal(theta[, 2], ga) > runif(n_item * n_person, 0, 1))
x <- 1 * (p_irt(theta[, 1], b, beta, w) > runif(n_item * n_person, 0, 1))

compute_normal_const <- function(ind, q_pts, q_w, x, w, b, beta, ga, sigma) {
  p_x <- p_irt(rep(q_pts[ind, 1], dim(w)[2]), b, beta, w)
  p_w <- p_cal(rep(q_pts[ind, 2], dim(w)[2]), ga)
  l_x <- x * log(p_x) + (1 - x) * log(1 - p_x)
  l_w <- w * log(p_w) + (1 - w) * log(1 - p_w)
  l <- apply(l_x, 2, sum) + apply(l_w, 2, sum) + q_w[ind]
  return(exp(l))
}

####################################################
############# M-Step ################
######## update ga ########
obj_ga <- function(ga, post_w_k, w_j, quad_pts) {
  out <- sum(post_w_k / (1 + exp(-outer(quad_pts, ga, "-")))) - w_j
  return(out)
}

update_ga <- function(i, post_w_k, w_j, quad_pts) {
  out <- uniroot(f = obj_ga, interval = c(-10, 10), post_w_k = post_w_k, w_j =
   w_j[i], quad_pts = quad_pts[, 2])$root
  return(out)
}

################# update rho ##############

obj_rho <- function(rho, expec_theta_squared, expec_theta_prod, n) {
  out <- n * rho * (1 - rho^2) + expec_theta_prod * (1 + rho^2) -
  expec_theta_squared * rho
  return(out)
}
############# update beta ################

obj_b <- function(b, quad_pts, post_w, x_j, w) {
  post_w <- post_w[w == 0, ]
  p_j_expec <- sum(p_irt_no_w(b, quad_pts) %*% t(post_w))
  out <- p_j_expec - x_j
  return(out)
}

update_b <- function(i, quad_pts, post_w, x_j, w) {
  out <- uniroot(obj_b, interval = c(-10, 10), quad_pts = quad_pts[, 1],
   post_w = post_w, x_j = x_j[i, 1], w = w[i, ])$root
  return(out)
}
###########################################

############## update beta ################
p_irt_w <- function(beta, theta, b) {
  out <- 1 / (1 + exp(outer(b, theta, "-") - beta))
  return(out)
}

obj_beta <- function(beta, b, quad_pts, post_w, x_j, w) {
  post_w <- post_w[w == 1, ]
  p_j_expec <- sum(p_irt_w(beta, quad_pts, b) %*% t(post_w))
  out <- p_j_expec - x_j
  return(out)
}

update_beta <- function(i, b, quad_pts, post_w, x_j, w) {
  out <- uniroot(obj_beta, interval = c(-10, 10), b = b[i], quad_pts =
    quad_pts[, 1], post_w = post_w, x_j = x_j[i, 2], w = w[i, ])$root
  return(out)
}
###########################################
########## true parameter values #########
true_b <- b
true_rho <- rho
true_ga <- ga
true_beta <- beta
###################################################################
########## compute some statistics and initialize parameters#######
x_j <- matrix(NA, nrow = n_item, ncol = 2)
x_j[, 2] <- apply(x * w, 1, sum) # sum score for w = 1
x_j[, 1] <- apply(x * (1 - w), 1, sum) # sum score for w = 0
w_j <- apply(w, 1, sum)
rho <- 0.5
b <- rnorm(n_item)
beta <- rnorm(n_item)
ga <- rnorm(n_item)
temp <- MultiGHQuad::init.quad(Q = 2, prior = list(mu = c(0, 0), Sigma = sig),
 ip = 25)
quad_w <- temp$W #log weights
quad_pts <- temp$X

############## EM iterations ####################
for (i in 1:25) {
  ###### compute posterior weights for each person ########
  post_w <- sapply(seq_len(length(quad_w)), FUN = compute_normal_const, q_pts
    = quad_pts, q_w = quad_w, x = x, w = w, b = b, beta = beta, ga = ga,
    sigma = sig)
  z <- apply(post_w, 1, sum)
  post_w <- post_w / z
  post_w_k <- apply(post_w, 2, sum)
  ga <- sapply(1:n_item, update_ga, post_w_k = post_w_k, w_j = w_j, quad_pts =
   quad_pts)
  b <- sapply(1:n_item, update_b, quad_pts = quad_pts, post_w = post_w, x_j =
   x_j, w = w)
  beta <- sapply(1:n_item, update_beta, b = b, quad_pts = quad_pts, post_w =
   post_w, x_j = x_j, w = w)
  expec_theta_squared <- sum(apply(post_w %*% (quad_pts^2), 2, sum))
  expec_theta_prod <- sum(post_w %*% apply(quad_pts, 1, prod))
  rho <- uniroot(f = obj_rho, interval = c(-1, 1), expec_theta_squared =
   expec_theta_squared, expec_theta_prod = expec_theta_prod, n = n_person)$root
  print(i)
}
