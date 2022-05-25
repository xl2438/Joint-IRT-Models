#include <RcppEigen.h>
#include <Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_poly.h>
#include <stdlib.h>
#include <stdio.h>
#include "quadrature.hpp"

// [[Rcpp::depends(RcppGSL)]]

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector2d;
using Eigen::Vector3d;


MatrixXd outer_substract(const VectorXd v1, const VectorXd v2) {
  MatrixXd p;
  p = v1.rowwise().replicate(v2.size()) - v2.transpose().colwise().replicate
  (v1.size());
  return p;
}

MatrixXd irt_2pl(VectorXd theta, VectorXd a, VectorXd b) {
  MatrixXd p;
	p = outer_substract(theta, b).cwiseProduct(-1 * a.transpose().colwise
    ().replicate(theta.size()));
  p = (1 + p.array().exp()).matrix().cwiseInverse();
  return p.transpose(); //J by N
}

MatrixXd p_irt(VectorXd theta, VectorXd a, VectorXd b,
  VectorXd ga, MatrixXd w) {
  MatrixXd p;
  p = (outer_substract(b, theta) - ga.rowwise().replicate(w.cols()).cwiseProduct
    (w)).cwiseProduct(a.rowwise().replicate(theta.size()));
  p = (1 + p.array().exp()).matrix().cwiseInverse();
  return p;
}

// [[Rcpp::export]]
MatrixXd compute_posterior(MatrixXd q_pts, VectorXd q_w,
  MatrixXd x, MatrixXd w, VectorXd a, VectorXd b,
  VectorXd alpha, VectorXd beta, VectorXd ga) {
  MatrixXd pw, mat_l, px10, px11, mat_w;
  VectorXd theta;
  theta = q_pts.col(1);
  pw = irt_2pl(theta, alpha, beta);

  mat_w.resize(a.size(), q_w.size());
  mat_w.fill(1);
  theta = q_pts.col(0);
  px11 = p_irt(theta, a, b, ga, mat_w).transpose();//QXJ  
  px10 = p_irt(theta, a, b, ga, (1-mat_w.array()).matrix()).transpose();

  mat_l = px11.array().log().matrix() * (x.cwiseProduct(w)) +
    (1-px10.array()).log().matrix() * (1-x.array().pow(1-w.array())).matrix() +
    (1-px11.array()).log().matrix() * (1 - x.array().pow(w.array())).matrix() +
    px10.array().log().matrix() * (1 - w.array().pow(x.array())).matrix();
  mat_l += q_w.rowwise().replicate(x.cols());
  mat_l += (pw.array().log().matrix().transpose() * w + (1-pw.array
      ()).log().matrix().transpose() * (1-w.array()).matrix());
  mat_l = mat_l.array().exp().matrix();
  mat_l = mat_l.cwiseProduct(mat_l.colwise().sum().colwise().replicate(q_w.size
    ()).cwiseInverse());
  return mat_l; //returns a Q X N matrix
}

gsl_vector grad_alpha_beta(VectorXd par, VectorXd w, VectorXd quad_theta,
 MatrixXd post_w) {
  VectorXd p, alpha(1), beta(1);
  MatrixXd res;
  gsl_vector * grad = gsl_vector_alloc(2);
  alpha(0) = par(0);
  beta(0) = par(1);
  p = irt_2pl(quad_theta, alpha, beta).transpose();
  res = outer_substract(w, p).cwiseProduct(post_w.transpose()); 
  gsl_vector_set(grad, 0, (res * (quad_theta.array() - beta(0)).matrix()).sum
    ());
  gsl_vector_set(grad, 1, -res.sum());
  return *grad;
}

gsl_matrix jac_alpha_beta(VectorXd par, VectorXd w, VectorXd quad_theta,
 MatrixXd post_w) {
  VectorXd p, alpha(1), beta(1), q, pq;
  MatrixXd res_weighted, diff_theta_beta, diff_theta_beta_pq;
  gsl_matrix * J = gsl_matrix_alloc(2,2);
  alpha(0) = par(0);
  beta(0) = par(1);
  p = irt_2pl(quad_theta, alpha, beta).transpose();
  q = (1 - p.array()).matrix();
  pq = p.cwiseProduct(q);
  res_weighted = outer_substract(w, p).cwiseProduct(post_w.transpose());
  diff_theta_beta = (quad_theta.array() - beta(0)).matrix();
  diff_theta_beta_pq = diff_theta_beta.cwiseProduct(pq);
  gsl_matrix_set(J, 0, 0, -(post_w.transpose() * 
    diff_theta_beta.cwiseProduct(diff_theta_beta_pq)).sum());
  gsl_matrix_set(J, 1, 0, (post_w.transpose() * diff_theta_beta_pq).sum());
  gsl_matrix_set(J, 0, 1, alpha(0) * gsl_matrix_get(J, 1, 0) -
    res_weighted.sum());
  gsl_matrix_set(J, 1, 1, -(alpha * pq.transpose() * post_w).sum());
  return *J;
}
 
struct rparams{
  VectorXd x;
  VectorXd w;
  VectorXd quad_theta;
  MatrixXd post_w;
};

int alpha_beta_f(const gsl_vector * x_par, void *params, gsl_vector * f) {
  VectorXd w = ((struct rparams *) params) -> w;
  VectorXd quad_theta = ((struct rparams *) params) -> quad_theta;
  MatrixXd post_w = ((struct rparams *) params) -> post_w;
  Vector2d par(gsl_vector_get(x_par, 0), gsl_vector_get(x_par, 1));
  gsl_vector temp = grad_alpha_beta(par, w, quad_theta, post_w);
  gsl_vector_memcpy(f, &temp);
  return GSL_SUCCESS;
}

int alpha_beta_df(const gsl_vector * x_par, void *params, gsl_matrix * J) {
  const VectorXd w = ((struct rparams *) params) -> w;
  const VectorXd quad_theta = ((struct rparams *) params) -> quad_theta;
  const MatrixXd post_w = ((struct rparams *) params) -> post_w;
  const Vector2d par(gsl_vector_get(x_par, 0), gsl_vector_get(x_par, 1));
  gsl_matrix temp = jac_alpha_beta(par, w, quad_theta, post_w);
  gsl_matrix_memcpy(J, &temp);
  return GSL_SUCCESS;
}

int alpha_beta_fdf(const gsl_vector * x_par, void *params, gsl_vector * f,
  gsl_matrix * J) {
  VectorXd w = ((struct rparams *) params) -> w;
  VectorXd quad_theta = ((struct rparams *) params) -> quad_theta;
  MatrixXd post_w = ((struct rparams *) params) -> post_w;

  VectorXd p, alpha(1), beta(1), q, pq;
  MatrixXd res_weighted, diff_theta_beta, diff_theta_beta_pq;
  double res_weighted_sum;
  alpha(0) = gsl_vector_get(x_par, 0);
  beta(0) = gsl_vector_get(x_par, 1);
  // compute involved terms //////
  p = irt_2pl(quad_theta, alpha, beta).transpose();
  q = (1 - p.array()).matrix();
  pq = p.cwiseProduct((1-p.array()).matrix());
  res_weighted = outer_substract(w, p).cwiseProduct(post_w.transpose());
  diff_theta_beta = (quad_theta.array() - beta(0)).matrix();
  res_weighted_sum = res_weighted.sum();
  diff_theta_beta_pq = diff_theta_beta.cwiseProduct(pq);
  ////////////////////////////////
  // compute the gradient vector//
  gsl_vector_set(f, 0, (res_weighted * diff_theta_beta).sum());
  gsl_vector_set(f, 1, -res_weighted_sum);
  /////////////////////////////////
  // compute the jacobian matrix ///
  gsl_matrix_set(J, 0, 0, -(post_w.transpose() * 
    diff_theta_beta.cwiseProduct(diff_theta_beta_pq)).sum());
  gsl_matrix_set(J, 1, 0, (post_w.transpose() * diff_theta_beta_pq).sum());
  gsl_matrix_set(J, 0, 1, alpha(0) * gsl_matrix_get(J, 1, 0) - res_weighted_sum);
  gsl_matrix_set(J, 1, 1, -(alpha * pq.transpose() * post_w).sum());
  ///////////////////////////////////
  return GSL_SUCCESS;
}

gsl_vector grad_a_b_gamma(VectorXd par, VectorXd x, VectorXd w, VectorXd
  quad_theta, MatrixXd post_w) {
  VectorXd p0, p1, a(2), b(2), gamma(2);
  MatrixXd res, p, theta_b_w_gamma, w_mat_const(2, quad_theta.size()), w_mat
  (w.size(), 2);
  gsl_vector * grad = gsl_vector_alloc(3);
  a.fill(par(0));
  b.fill(par(1));
  gamma.fill(par(2));
  w_mat_const.row(0).fill(0);
  w_mat_const.row(1).fill(1);
  w_mat.col(0) = (1-w.array()).matrix();
  w_mat.col(1) = w;
  p =  (w_mat * p_irt(quad_theta, a, b, gamma, w_mat_const)).transpose();
  res = (x.transpose().colwise().replicate(quad_theta.size()) - p).cwiseProduct
  (post_w); //Q X N
  theta_b_w_gamma = outer_substract((quad_theta.array()-b(0)).matrix(),
    -gamma(0)*w); //Q X N
  gsl_vector_set(grad, 0, res.cwiseProduct(theta_b_w_gamma).sum());
  gsl_vector_set(grad, 1, -res.sum());
  gsl_vector_set(grad, 2, (res*w).sum());
  return *grad;
}

gsl_matrix jac_a_b_gamma(VectorXd par, VectorXd x, VectorXd w, VectorXd
  quad_theta, MatrixXd post_w) {
  VectorXd p0, p1, a(2), b(2), gamma(2);
  MatrixXd res, p, q, pq, theta_b_w_gamma, w_mat_const(2, quad_theta.size()),
  w_mat(w.size(), 2),  theta_b_w_gamma_pq_weighted;
  gsl_matrix * J = gsl_matrix_alloc(3, 3);
  a.fill(par(0));
  b.fill(par(1));
  gamma.fill(par(2));
  w_mat_const.row(0).fill(0);
  w_mat_const.row(1).fill(1);
  w_mat.col(0) = (1-w.array()).matrix();
  w_mat.col(1) = w;
  p =  (w_mat * p_irt(quad_theta, a, b, gamma, w_mat_const)).transpose();
  q = (1 - p.array()).matrix();
  pq = p.cwiseProduct(q).cwiseProduct(post_w);
  res = (x.transpose().colwise().replicate(quad_theta.size()) - p).cwiseProduct
  (post_w); //Q X N
  theta_b_w_gamma = outer_substract((quad_theta.array()-b(0)).matrix(),
    -gamma(0)*w); //Q X N
  theta_b_w_gamma_pq_weighted = theta_b_w_gamma.cwiseProduct(pq);
  double theta_b_w_gamma_pq_weighted_sum = theta_b_w_gamma_pq_weighted.sum();
  
  gsl_matrix_set(J, 0, 0, -theta_b_w_gamma_pq_weighted.cwiseProduct
    (theta_b_w_gamma).sum());
  gsl_matrix_set(J, 0, 1, a(0)*theta_b_w_gamma_pq_weighted_sum - res.sum());
  gsl_matrix_set(J, 0, 2, -(a(0)*theta_b_w_gamma_pq_weighted*w).sum() + 
    (res*w).sum());
  gsl_matrix_set(J, 1, 0, theta_b_w_gamma_pq_weighted_sum);
  gsl_matrix_set(J, 1, 1, -a(0)*pq.sum());
  gsl_matrix_set(J, 1, 2, a(0)*(pq*w).sum());
  gsl_matrix_set(J, 2, 0, -(theta_b_w_gamma_pq_weighted*w).sum());
  gsl_matrix_set(J, 2, 1, gsl_matrix_get(J, 1, 2));
  gsl_matrix_set(J, 2, 2, -gsl_matrix_get(J, 1, 2));

  return *J;
}

int a_b_gamma_f(const gsl_vector * x_par, void *params, gsl_vector * f) {
  VectorXd x = ((struct rparams *) params) -> x;
  VectorXd w = ((struct rparams *) params) -> w;
  VectorXd quad_theta = ((struct rparams *) params) -> quad_theta;
  MatrixXd post_w = ((struct rparams *) params) -> post_w;
  Vector3d par(gsl_vector_get(x_par, 0), gsl_vector_get(x_par, 1),
   gsl_vector_get(x_par, 2));
  gsl_vector temp = grad_a_b_gamma(par, x, w, quad_theta, post_w);
  gsl_vector_memcpy(f, &temp);
  return GSL_SUCCESS;
}

int a_b_gamma_df(const gsl_vector * x_par, void *params, gsl_matrix * J) {
  VectorXd x = ((struct rparams *) params) -> x;
  VectorXd w = ((struct rparams *) params) -> w;
  VectorXd quad_theta = ((struct rparams *) params) -> quad_theta;
  MatrixXd post_w = ((struct rparams *) params) -> post_w;
  Vector3d par(gsl_vector_get(x_par, 0), gsl_vector_get(x_par, 1),
   gsl_vector_get(x_par, 2));
  gsl_matrix temp = jac_a_b_gamma(par, x, w, quad_theta, post_w);
  gsl_matrix_memcpy(J, &temp);
  return GSL_SUCCESS;  
}

int a_b_gamma_fdf(const gsl_vector * x_par, void *params, gsl_vector * f,
  gsl_matrix * J) {
  VectorXd x = ((struct rparams *) params) -> x;
  VectorXd w = ((struct rparams *) params) -> w;
  VectorXd quad_theta = ((struct rparams *) params) -> quad_theta;
  MatrixXd post_w = ((struct rparams *) params) -> post_w;

  VectorXd a(1), b(1), gamma(1);
  MatrixXd w_mat_const(2, quad_theta.size()), w_mat(w.size(), 2),
  theta_b_w_gamma_pq_weighted, res, p, q, pq, theta_b_w_gamma;
  double theta_b_w_gamma_pq_weighted_sum, res_sum,
  theta_b_w_gamma_pq_weighted_w_sum, res_w_sum;

  a.fill(gsl_vector_get(x_par, 0));
  b.fill(gsl_vector_get(x_par, 1));
  gamma.fill(gsl_vector_get(x_par, 2));
  w_mat_const.row(0).fill(0);
  w_mat_const.row(1).fill(1);
  w_mat.col(0) = (1-w.array()).matrix();
  w_mat.col(1) = w;
  p =  (w_mat * p_irt(quad_theta, a, b, gamma, w_mat_const)).transpose();
  q = (1 - p.array()).matrix();
  pq = p.cwiseProduct(q).cwiseProduct(post_w);
  res = (x.transpose().colwise().replicate(quad_theta.size()) - p).cwiseProduct
  (post_w); //Q X N
  theta_b_w_gamma = outer_substract((quad_theta.array()-b(0)).matrix(),
    -gamma(0)*w); //Q X N
  theta_b_w_gamma_pq_weighted = theta_b_w_gamma.cwiseProduct(pq);
  theta_b_w_gamma_pq_weighted_sum = theta_b_w_gamma_pq_weighted.sum();
  res_sum = res.sum();
  theta_b_w_gamma_pq_weighted_w_sum = (theta_b_w_gamma_pq_weighted*w).sum();
  res_w_sum = (res*w).sum();

  gsl_vector_set(f, 0, res.cwiseProduct(theta_b_w_gamma).sum());
  gsl_vector_set(f, 1, -res_sum);
  gsl_vector_set(f, 2, res_w_sum);

  gsl_matrix_set(J, 0, 0, -theta_b_w_gamma_pq_weighted.cwiseProduct
    (theta_b_w_gamma).sum());
  gsl_matrix_set(J, 0, 1, a(0)*theta_b_w_gamma_pq_weighted_sum - res_sum);
  gsl_matrix_set(J, 0, 2, -a(0)*theta_b_w_gamma_pq_weighted_w_sum + 
    res_w_sum);
  gsl_matrix_set(J, 1, 0, theta_b_w_gamma_pq_weighted_sum);
  gsl_matrix_set(J, 1, 1, -a(0)*pq.sum());
  gsl_matrix_set(J, 1, 2, a(0)*(pq*w).sum());
  gsl_matrix_set(J, 2, 0, -theta_b_w_gamma_pq_weighted_w_sum);
  gsl_matrix_set(J, 2, 1, gsl_matrix_get(J, 1, 2));
  gsl_matrix_set(J, 2, 2, -gsl_matrix_get(J, 1, 2));

  return GSL_SUCCESS;
}

struct rparams_mat
{
  MatrixXd x;
  MatrixXd w;
  MatrixXd quad_theta;
  MatrixXd post_w;
};

void solve(gsl_multiroot_fdfsolver *s) {
  int status;
  size_t iter = 0;
  do{
    iter++;
    status = gsl_multiroot_fdfsolver_iterate(s);
    if(status) break;
    status = gsl_multiroot_test_delta(s->dx, s->x, 1e-6, 1e-7);
  }while(status == GSL_CONTINUE && iter < 150);
}

void estimate_params(gsl_matrix * par_a_b_gamma, gsl_matrix * par_alpha_beta,
  double* rho, void *params) {
  MatrixXd x = ((struct rparams_mat *) params) -> x;
  MatrixXd w = ((struct rparams_mat *) params) -> w;
  MatrixXd quad_theta = ((struct rparams_mat *) params) -> quad_theta;
  MatrixXd post_w = ((struct rparams_mat *) params) -> post_w;

  const gsl_multiroot_fdfsolver_type *T;
  gsl_multiroot_fdfsolver *s;
  gsl_multiroot_function_fdf f;
  T = gsl_multiroot_fdfsolver_hybridsj;
  for (int i = 0; i < x.rows(); i++){ 
    gsl_vector * x_par = gsl_vector_alloc(3);
    gsl_matrix_get_row(x_par, par_a_b_gamma, i);
    struct rparams p = {x.row(i), w.row(i), quad_theta.col(0), post_w};
    f = {&a_b_gamma_f, &a_b_gamma_df, &a_b_gamma_fdf, 3, &p};
    s = gsl_multiroot_fdfsolver_alloc(T, 3);
    gsl_multiroot_fdfsolver_set(s, &f, x_par);
    solve(s);
    gsl_matrix_set_row(par_a_b_gamma, i, s->x);
    gsl_vector_free(x_par);
    gsl_multiroot_fdfsolver_free(s);

    p = {x.row(i), w.row(i), quad_theta.col(1), post_w};
    x_par = gsl_vector_alloc(2);
    gsl_matrix_get_row(x_par, par_alpha_beta, i);
    f = {&alpha_beta_f, &alpha_beta_df, &alpha_beta_fdf, 2, &p};
    s = gsl_multiroot_fdfsolver_alloc(T, 2);
    gsl_multiroot_fdfsolver_set(s, &f, x_par);
    solve(s);
    gsl_matrix_set_row(par_alpha_beta, i, s->x);
    gsl_vector_free(x_par);
    gsl_multiroot_fdfsolver_free(s);
  }

  double a, b, c;
  a = -(quad_theta.col(0).cwiseProduct(quad_theta.col(1)).transpose() *
    post_w).sum() / post_w.cols();
  b = ((quad_theta.col(0).cwiseProduct(quad_theta.col(0)) + quad_theta.col
    (1).cwiseProduct(quad_theta.col(1))).transpose() * post_w).sum()
  / post_w.cols() - 1;
  c = a;
  double* rho0 = new double;
  double* rho1 = new double;
  double* rho2 = new double;
  int n_root;
  n_root = gsl_poly_solve_cubic(a, b, c, rho0, rho1, rho2);
  if (n_root == 1) *rho = *rho0;
  else if (1 <= *rho0 && *rho0 >= -1) *rho = *rho0;
  else if (1 <= *rho1 && *rho1 >= -1) *rho = *rho1;
  else
    *rho = *rho2;
  delete rho0;
  delete rho1;
  delete rho2;
}

// [[Rcpp::export]]
double update_a_b_gamma(MatrixXd x, MatrixXd w, MatrixXd quad_theta, 
  MatrixXd post_w) {
  gsl_matrix * par_a_b_gamma = gsl_matrix_alloc(x.rows(), 3);
  gsl_matrix * par_alpha_beta = gsl_matrix_alloc(x.rows(), 2);
  for (int i = 0; i < x.rows(); i++) {
    gsl_matrix_set(par_a_b_gamma, i, 0, 1.15);
    gsl_matrix_set(par_a_b_gamma, i, 1, 0.50);
    gsl_matrix_set(par_a_b_gamma, i, 2, 0.45);
    gsl_matrix_set(par_alpha_beta, i, 0, 1.15);
    gsl_matrix_set(par_alpha_beta, i, 1, 0.50);
  }
  struct rparams_mat p = {x, w, quad_theta, post_w};
  double* rho = new double;
  *rho = 2.0;
  estimate_params(par_a_b_gamma, par_alpha_beta, rho, &p);
  return *rho;
}

// [[Rcpp::export]]
MatrixXd compute_quad_pts(double rho, int ip) {
  Quadrature my_quadrature(rho, ip);
  return my_quadrature.get_quad_pts();
} 

