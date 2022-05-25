#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector2d;

class quadrature {
private:
  MatrixXd quad_pts;
  VectorXd quad_w; //log weights
  double rho;
  int ip;

public:
  quadrature(int ip = 13) {
    this->ip = ip;
    rho = 0.2;
  }

  void update(double rho) {
    this->rho = rho;
    MatrixXd sigma(2,2);
    sigma(0, 0) = 1.0;
    sigma(0, 1) = rho;
    sigma(1, 0) = rho;
    sigma(1, 1) = 1.0;
    Rcpp::Environment pkg = Rcpp::Environment::namespace_env("MultiGHQuad");
    Rcpp::Function f = pkg["init.quad"];
    Vector2d mu(0.0, 0.0);
    Rcpp::List prior = Rcpp::List::create(Rcpp::Named("mu") = mu, Rcpp::Named
      ("Sigma") = sigma);
    Rcpp::List result_f = f(Rcpp::Named("prior", prior), Rcpp::Named("ip", ip));
    quad_pts = result_f[0];
    quad_w = result_f[1];
  }

  MatrixXd get_quad_pts(){return quad_pts;}
  VectorXd get_quad_w(){return quad_w;}
};

#endif