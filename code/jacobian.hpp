#ifndef JACOBIAN_HPP
#define	JACOBIAN_HPP

#include "quadrature.hpp"
#include "jointModel.hpp"
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::vectorXd;

MatrixXd hessian_complete_data_loglik_alpha_beta(jointModel& mo, VectorXd theta) {
  VectorXd alpha, beta;
  MatrixXd w;
  alpha = mo.getAlpha();
  beta = mo.getBeta();
  w = mo.getW();
  
}

#endif