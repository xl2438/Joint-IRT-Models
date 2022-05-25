#ifndef DERIVATIVES_HPP
#define DERIVATIVES_HPP

#include <RcppGSL.h>
#include <RcppEigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "utilities.hpp"

// [[Rcpp::depends(RcppGSL)]]

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector2d;
using Eigen::Vector3d;

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

  VectorXd p, alpha(1), beta(1);
  MatrixXd res;
  alpha(0) = par(0);
  beta(0) = par(1);
  p = irt_2pl(quad_theta, alpha, beta).transpose();
  res = outer_substract(w, p).cwiseProduct(post_w.transpose()); 
  gsl_vector_set(f, 0, (res * (quad_theta.array() - beta(0)).matrix()).sum
    ());
  gsl_vector_set(f, 1, -res.sum());
  return GSL_SUCCESS;
}

int alpha_beta_df(const gsl_vector * x_par, void *params, gsl_matrix * J) {
  VectorXd w = ((struct rparams *) params) -> w;
  VectorXd quad_theta = ((struct rparams *) params) -> quad_theta;
  MatrixXd post_w = ((struct rparams *) params) -> post_w;
  Vector2d par(gsl_vector_get(x_par, 0), gsl_vector_get(x_par, 1));

  VectorXd p, alpha(1), beta(1), q, pq;
  MatrixXd res_weighted, diff_theta_beta, diff_theta_beta_pq;
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

int a_b_gamma_f(const gsl_vector * x_par, void *params, gsl_vector * f) {
  VectorXd x = ((struct rparams *) params) -> x;
  VectorXd w = ((struct rparams *) params) -> w;
  VectorXd quad_theta = ((struct rparams *) params) -> quad_theta;
  MatrixXd post_w = ((struct rparams *) params) -> post_w;
  Vector3d par(gsl_vector_get(x_par, 0), gsl_vector_get(x_par, 1),
   gsl_vector_get(x_par, 2));
  
  VectorXd p0, p1, a(2), b(2), gamma(2);
  MatrixXd res, p, theta_b_w_gamma, w_mat_const(2, quad_theta.size()), w_mat
  (w.size(), 2);
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
  gsl_vector_set(f, 0, res.cwiseProduct(theta_b_w_gamma).sum());
  gsl_vector_set(f, 1, -res.sum());
  gsl_vector_set(f, 2, (res*w).sum());
  return GSL_SUCCESS;
}

int a_b_gamma_df(const gsl_vector * x_par, void *params, gsl_matrix * J) {
  VectorXd x = ((struct rparams *) params) -> x;
  VectorXd w = ((struct rparams *) params) -> w;
  VectorXd quad_theta = ((struct rparams *) params) -> quad_theta;
  MatrixXd post_w = ((struct rparams *) params) -> post_w;
  Vector3d par(gsl_vector_get(x_par, 0), gsl_vector_get(x_par, 1),
   gsl_vector_get(x_par, 2));
  
  VectorXd p0, p1, a(2), b(2), gamma(2);
  MatrixXd res, p, q, pq, theta_b_w_gamma, w_mat_const(2, quad_theta.size()),
  w_mat(w.size(), 2),  theta_b_w_gamma_pq_weighted;
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

#endif