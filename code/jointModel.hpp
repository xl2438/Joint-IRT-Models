#ifndef JOINTMODEL_HPP
#define JOINTMODEL_HPP

#include "model.hpp"
#include <RcppEigen.h>
#include "utilities_rand.hpp"

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;

struct T_data {
  MatrixXd x;
  MatrixXd w;
};

struct T_par {
  VectorXd a;
  VectorXd b;
  VectorXd ga;
  VectorXd alpha;
  VectorXd beta;
  double rho;
};

class jointModel: public model<T_data, T_par> {
  public:
    jointModel(MatrixXd& x, MatrixXd& w);
    int updatePar(T_par new_par);
    T_par getPar(){return par;};
    T_data getData(){return data;};
    MatrixXd getX(){return data.x;};
    MatrixXd getW(){return data.w;};
    VectorXd getA(){return par.a;};
    VectorXd getB(){return par.b;};
    VectorXd getGa(){return par.ga;};
    VectorXd getAlpha(){return par.alpha;};
    VectorXd getBeta(){return par.beta;};
    double getRho(){return par.rho;};
};

jointModel::jointModel(MatrixXd& x, MatrixXd& w) {
  data.x = x;
  data.w = w;
  par.a.resize(data.x.rows());
  par.b.resize(data.x.rows());
  par.a.fill(rlnorm(0, 0.15));
  par.b.fill(rnorm(0, 1));
  par.alpha = par.a;
  par.beta = par.b;
  par.ga.resize(data.x.rows());
  par.ga.fill(rnorm(0, 1));
  par.rho = 0.2;
}

int jointModel::updatePar(T_par new_par) {
  par.a = new_par.a;
  par.b = new_par.b;
  par.ga = new_par.ga;
  par.alpha = new_par.alpha;
  par.beta = new_par.beta;
  par.rho = new_par.rho;
  return 0;
}

#endif