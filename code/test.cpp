#include <RcppEigen.h>
#include <Rcpp.h>
#include "jointModel.hpp"
#include "em.hpp"
#include "quadrature.hpp"
#include "posterior.hpp"
#include "se.hpp"
#include <iostream>
#include <string>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;

// [[Rcpp::export]]
Rcpp::List estimate(MatrixXd& x, MatrixXd& w, int n_itr = 80) {
  jointModel mo(x, w);
  quadrature quad(13);
  MatrixXd post_w;

  int i = 0;
  double previous_loglik = 0;
  double current_loglik = 0;
  do{
  	previous_loglik = current_loglik;
  	quad.update(mo.getRho());
  	post_w = compute_posterior(mo, quad);
  	updateParams(mo, quad, post_w);
  	current_loglik = updateLoglik(mo, quad);
  	std::cout << "Iteration " << i+1 << " Log-likelihood:" << current_loglik
  	<< std::endl;
  	i++;
  } while(i < n_itr && abs(current_loglik - previous_loglik) > 1e-1);

  
  return Rcpp::List::create(Rcpp::Named("a") = mo.getA(), Rcpp::Named("b") =
  	mo.getB(), Rcpp::Named("alpha") = mo.getAlpha(), Rcpp::Named("beta") =
  	mo.getBeta(), Rcpp::Named("ga") = mo.getGa(), Rcpp::Named("rho") = mo.getRho()
  	// , Rcpp::Named("se") = computeSE(mo)
    // , Rcpp::Named("x") = x, Rcpp::Named("w") = w
  	);
}