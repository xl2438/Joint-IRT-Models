#ifndef SE_HPP
#define SE_HPP

#include <RcppEigen.h>
#include <Rcpp.h>
#include "quadrature.hpp"
#include "jointModel.hpp"

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector2d;

// [[Rcpp::export]]
double obsLoglik(VectorXd& params, MatrixXd& dat_x, MatrixXd& w, int j) {
  VectorXd alpha = params.segment(0, j);
  VectorXd beta = params.segment(j, j);
  VectorXd a = params.segment(2*j, j);
  VectorXd b = params.segment(3*j, j);
  VectorXd ga = params.segment(4*j, j);
  double rho = params.segment(5*j, 1)(0);
  MatrixXd x = dat_x;

  quadrature quad;
  quad.update(rho);
  MatrixXd q_pts = quad.get_quad_pts();
  MatrixXd q_w = quad.get_quad_w();

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

MatrixXd computeLoglik(jointModel& mo) {
  MatrixXd x = mo.getX();
  MatrixXd w = mo.getW();
  VectorXd a = mo.getA();
  VectorXd b = mo.getB();
  VectorXd ga = mo.getGa();
  VectorXd alpha = mo.getAlpha();
  VectorXd beta = mo.getBeta();
  double rho = mo.getRho();

  VectorXd params(5*alpha.size()+1);
  params << alpha, beta, a, b, ga, rho;

  Rcpp::Function f("obsLoglik");
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("numDeriv");
  Rcpp::Function g = pkg["hessian"];

  Rcpp::NumericMatrix out =  g(Rcpp::Named("func", f), Rcpp::Named("x",
   params), Rcpp::Named("dat_x",x), Rcpp::Named("w",w), Rcpp::Named("j", a.size
   ()));

  return Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(out);
}

VectorXd computeSE(jointModel& mo) {
  MatrixXd info = computeLoglik(mo);
  info = -info;
  VectorXd se;
  se = info.inverse().diagonal().cwiseSqrt();
  return se;
}

#endif