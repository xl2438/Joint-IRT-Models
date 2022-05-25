#ifndef EM_HPP
#define EM_HPP

#include <RcppEigen.h>
#include <Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_poly.h>
#include "quadrature.hpp"
#include "derivatives.hpp"
#include "jointModel.hpp"
#include "utilities_rand.hpp"
#include <iostream>

// [[Rcpp::depends(RcppGSL)]]

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector2d;
using Eigen::Vector3d;

int solve(gsl_multiroot_fdfsolver *s) {
  int status;
  size_t iter = 0;
  do{
    iter++;
    status = gsl_multiroot_fdfsolver_iterate(s);
    if(status) break;
    // status = gsl_multiroot_test_delta(s->dx, s->x, 1e-11, 0);
    status = gsl_multiroot_test_residual(s->f, 1e-3);
  }while(status == GSL_CONTINUE && iter < 50);
  if(iter >= 50) std::cout << "Maximum Number of Iterations Reached!" <<
   std::endl;
  // std::cout << "# of iterations: " << iter << std::endl;
  return status;
}

int updateParams(jointModel& mo, quadrature& quad, MatrixXd& post_w) {
  MatrixXd x = mo.getData().x;
  MatrixXd w = mo.getData().w;
  MatrixXd quad_theta = quad.get_quad_pts();
  T_par par = mo.getPar();
  int status; 
  const gsl_multiroot_fdfsolver_type *T;
  gsl_multiroot_fdfsolver *s;
  gsl_multiroot_function_fdf f;
  T = gsl_multiroot_fdfsolver_hybridj;
  const gsl_multiroot_fdfsolver_type *G;
  G = gsl_multiroot_fdfsolver_gnewton;
  for (int i = 0; i < x.rows(); i++){ 
    gsl_vector * x_par = gsl_vector_alloc(3);
    gsl_vector_set(x_par, 0, par.a(i));
    gsl_vector_set(x_par, 1, par.b(i));
    gsl_vector_set(x_par, 2, par.ga(i));
    struct rparams p = {x.row(i), w.row(i), quad_theta.col(0), post_w};
    f = {&a_b_gamma_f, &a_b_gamma_df, &a_b_gamma_fdf, 3, &p};
    s = gsl_multiroot_fdfsolver_alloc(T, 3);
    gsl_multiroot_fdfsolver_set(s, &f, x_par);
    int counter = 0;
    while(solve(s) != GSL_SUCCESS) {
      // if (counter == 1) {
      //   s = gsl_multiroot_fdfsolver_alloc(G, 3);
      // }
      gsl_vector_set(x_par, 0, rlnorm(0, 0.15));
      gsl_vector_set(x_par, 1, rnorm(0, 1));
      gsl_vector_set(x_par, 2, rnorm(0, 1));
      gsl_multiroot_fdfsolver_set(s, &f, x_par);
      status = solve(s);
      std::cout << counter + 2 << " attempt item " << i << " status: " << status
      << std::endl;
      counter++;
    }
    par.a(i) = gsl_vector_get(s->x, 0);
    par.b(i) = gsl_vector_get(s->x, 1);
    par.ga(i) = gsl_vector_get(s->x, 2);
    gsl_vector_free(x_par);
    gsl_multiroot_fdfsolver_free(s);

    p = {x.row(i), w.row(i), quad_theta.col(1), post_w};
    x_par = gsl_vector_alloc(2);
    gsl_vector_set(x_par, 0, par.alpha(i));
    gsl_vector_set(x_par, 1, par.beta(i));
    // gsl_vector_set(x_par, 0, 1.01);
    // gsl_vector_set(x_par, 1, 0.1);
    f = {&alpha_beta_f, &alpha_beta_df, &alpha_beta_fdf, 2, &p};
    s = gsl_multiroot_fdfsolver_alloc(T, 2);
    gsl_multiroot_fdfsolver_set(s, &f, x_par);
    status = solve(s);
    // std::cout << "item " << i << " status: " << status << std::endl;
    counter = 0;
    if(status != GSL_SUCCESS) {
      if (counter == 1) {
        s = gsl_multiroot_fdfsolver_alloc(G, 3);
      }
      gsl_vector_set(x_par, 0, rlnorm(0, 0.15));
      gsl_vector_set(x_par, 1, rnorm(0, 1));
      gsl_multiroot_fdfsolver_set(s, &f, x_par);
      status = solve(s);
      std::cout << "2nd attempt item " << i << " status: " << status <<
      std::endl;
      counter++;
    }
    par.alpha(i) = gsl_vector_get(s->x, 0);
    par.beta(i) = gsl_vector_get(s->x, 1);
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
  if (n_root == 1) par.rho = *rho0;
  else if (1 <= *rho0 && *rho0 >= -1) par.rho = *rho0;
  else if (1 <= *rho1 && *rho1 >= -1) par.rho = *rho1;
  else
    par.rho = *rho2;
  delete rho0;
  delete rho1;
  delete rho2;

  mo.updatePar(par);
  return 0;
}

double updateLoglik(jointModel& mo, quadrature& quad) {
  MatrixXd q_pts = quad.get_quad_pts();
  MatrixXd q_w = quad.get_quad_w();
  MatrixXd x = mo.getX();
  MatrixXd w = mo.getW();
  VectorXd a = mo.getA();
  VectorXd b = mo.getB();
  VectorXd ga = mo.getGa();
  VectorXd alpha = mo.getAlpha();
  VectorXd beta = mo.getBeta();

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

  return mat_l.colwise().sum().array().log().sum();
}

#endif