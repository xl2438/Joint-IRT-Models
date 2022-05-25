#ifndef POSTERIOR_HPP
#define POSTERIOR_HPP

#include <RcppEigen.h>
#include "utilities.hpp"
#include "quadrature.hpp"
#include "jointModel.hpp"

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;

MatrixXd compute_posterior(jointModel& mo, quadrature& quad) {

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
  mat_l = mat_l.cwiseProduct(mat_l.colwise().sum().colwise().replicate(q_w.size
    ()).cwiseInverse());
  return mat_l; //returns a Q X N matrix
}

#endif