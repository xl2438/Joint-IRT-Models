#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;

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

#endif