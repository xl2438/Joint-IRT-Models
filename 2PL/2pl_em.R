irt_2pl <- function(theta, a, b) {
  out <- 1 / (1 + exp(a * outer(b, theta, "-")))
  return(out)
}

n_person <- 2000
n_item <- 20
a_true <- rlnorm(n_item, 0, 0.15)
b_true <- rnorm(n_item)
theta <- rnorm(n_person)

p <- irt_2pl(theta, a_true, b_true)
x <- 1 * (p > runif(length(p), 0, 1))

compute_joint <- function(ind, q_pts, x, a, b, log_w_pts) {
  p <- as.vector(irt_2pl(q_pts[ind], a, b))
  l <- x * log(p) + (1 - x) * log(1 - p)
  l <- apply(l, 2, sum) + log_w_pts[ind]
  return(exp(l))
}
compute_post_w <- function(q_pts, x, log_w_pts, a, b) {
  l <- sapply(seq_len(length(q_pts)), FUN = compute_joint, q_pts =
    q_pts, x = x, a = a, b = b, log_w_pts = log_w_pts)
  post_w <- l / apply(l, 1, sum)
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
  sln <- nleqslv::nleqslv(x = c(1, 0), fn = obj_fun, jac = jacobian_obj, x_vec =
   x[ind, ], q_pts = q_pts, post_w = post_w, method = "Newton")$x
  return(sln)
}
estimate_item_pars <- function(x, q_pts, post_w) {
  out <- sapply(c(1:dim(x)[1]), FUN = estimate_par_j, x = x, q_pts = q_pts,
    post_w = post_w)
  return(out)
}
####################### initialize ########################
a <- rlnorm(n_item, 0, 0.15)
b <- rnorm(n_item)
q_pts <- seq(from = -4, to = 4, length.out = 25)
log_w_pts <- dnorm(q_pts, log = TRUE)
#######################################################
start <- Sys.time()
for (i in 1:25) {
  post_w <- compute_post_w(q_pts, x, log_w_pts, a, b)
  item_par <- estimate_item_pars(x, q_pts, post_w)
  a <- item_par[1, ]
  b <- item_par[2, ]
}
Sys.time() - start

###################### check derivatives ###################################
# start <- Sys.time()
# for (i in 1:1000) {
#   temp <- jacobian_obj(c(1, 0), x[1, ], q_pts, post_w)
# }
# Sys.time() - start

# start <- Sys.time()
# for (i in 1:1000) {
#   temp <- numDeriv::jacobian(func = obj_fun_j, c(1, 0), x_vec = x[1, ], q_pts =
#     q_pts, post_w = post_w)
# }
# Sys.time() - start

# jacobian_obj(c(1, 0), x[1, ], q_pts, post_w)
# numDeriv::jacobian(func = obj_fun_j, c(1, 0), x_vec = x[1, ], q_pts =
#     q_pts, post_w = post_w)
############################################################################