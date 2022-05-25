#include <RcppEigen.h>
#include <Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_poly.h>

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

struct rparams_mat
{
  MatrixXd x;
  MatrixXd w;
  MatrixXd quad_theta;
  MatrixXd post_w;
};

class Model {
private:
  //data
  MatrixXd x;
  MatrixXd w;
  //model parameters
  gsl_matrix * par_a_b_gamma;
  gsl_matrix * par_alpha_beta;
  double* rho;

  //private member functions
  MatrixXd outer_substract(VectorXd v1, VectorXd v2);
  MatrixXd irt_2pl(VectorXd theta, VectorXd a, VectorXd b);
  MatrixXd p_irt(VectorXd theta, VectorXd a, VectorXd b, VectorXd ga, MatrixXd
   w);

  void solve(gsl_multiroot_fdfsolver *s);

public:
  //Constructor
  Model(MatrixXd x, MatrixXd w);

  //public member functions
  MatrixXd compute_posterior(MatrixXd q_pts, VectorXd q_w);
  int a_b_gamma_f(const gsl_vector * x_par, void *params, gsl_vector * f);
  int a_b_gamma_df(const gsl_vector * x_par, void *params, gsl_matrix * J);
  int a_b_gamma_fdf(const gsl_vector * x_par, void *params, gsl_vector * f,
    gsl_matrix * J);

  int alpha_beta_f(const gsl_vector * x_par, void *params, gsl_vector * f);
  int alpha_beta_df(const gsl_vector * x_par, void *params, gsl_matrix * J);
  int alpha_beta_fdf(const gsl_vector * x_par, void *params, gsl_vector * f,
    gsl_matrix * J);

  gsl_vector grad_a_b_gamma(VectorXd par, VectorXd x, VectorXd w, VectorXd
    quad_theta, MatrixXd post_w);
  gsl_matrix jac_a_b_gamma(VectorXd par, VectorXd x, VectorXd w, VectorXd
    quad_theta, MatrixXd post_w);
  gsl_vector grad_alpha_beta(VectorXd par, VectorXd w, VectorXd quad_theta,
    MatrixXd post_w);
  gsl_matrix jac_alpha_beta(VectorXd par, VectorXd w, VectorXd quad_theta,
    MatrixXd post_w);







  



};

Model::Model(MatrixXd x, MatrixXd w) {
  this->x = x;
  this->w = w;
  par_a_b_gamma = gsl_matrix_alloc(x.rows(), 3);
  par_alpha_beta = gsl_matrix_alloc(x.rows(), 2);
  for (int i = 0; i < x.rows(); i++) {
    gsl_matrix_set(par_a_b_gamma, i, 0, 1.15);
    gsl_matrix_set(par_a_b_gamma, i, 1, 0.50);
    gsl_matrix_set(par_a_b_gamma, i, 2, 0.45);
    gsl_matrix_set(par_alpha_beta, i, 0, 1.15);
    gsl_matrix_set(par_alpha_beta, i, 1, 0.50);
  }
  rho = new double;
  *rho = 0.2;
}

MatrixXd Model::outer_substract(VectorXd v1, VectorXd v2) {
  MatrixXd p;
  p = v1.rowwise().replicate(v2.size()) - v2.transpose().colwise().replicate
  (v1.size());
  return p;
}

MatrixXd Model::irt_2pl(VectorXd theta, VectorXd a, VectorXd b) {
  MatrixXd p;
  p = outer_substract(theta, b).cwiseProduct(-1 * a.transpose().colwise
    ().replicate(theta.size()));
  p = (1 + p.array().exp()).matrix().cwiseInverse();
  return p.transpose(); //J by N
}

MatrixXd Model::p_irt(VectorXd theta, VectorXd a, VectorXd b,
  VectorXd ga, MatrixXd w) {
  MatrixXd p;
  p = (outer_substract(b, theta) - ga.rowwise().replicate(w.cols()).cwiseProduct
    (w)).cwiseProduct(a.rowwise().replicate(theta.size()));
  p = (1 + p.array().exp()).matrix().cwiseInverse();
  return p;
}

MatrixXd Model::compute_posterior(MatrixXd q_pts, VectorXd q_w) {
  VectorXd a(x.rows()), b(x.rows()), ga(x.rows()), alpha(x.rows()), beta(x.rows
    ());
  for (int i = 0; i < x.rows(); i++) {
    a(i) = gsl_matrix_get(par_a_b_gamma, i, 0);
    b(i) = gsl_matrix_get(par_a_b_gamma, i, 1);
    ga(i) = gsl_matrix_get(par_a_b_gamma, i, 2);
    alpha(i) = gsl_matrix_get(par_alpha_beta, i, 0);
    beta(i) = gsl_matrix_get(par_alpha_beta, i, 1);
  }

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

//definitions of derivatives
gsl_vector Model::grad_alpha_beta(VectorXd par, VectorXd w, VectorXd quad_theta,
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

gsl_matrix Model::jac_alpha_beta(VectorXd par, VectorXd w, VectorXd quad_theta,
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
 
int Model::alpha_beta_f(const gsl_vector * x_par, void *params, gsl_vector * f)
{
  VectorXd w = ((struct rparams *) params) -> w;
  VectorXd quad_theta = ((struct rparams *) params) -> quad_theta;
  MatrixXd post_w = ((struct rparams *) params) -> post_w;
  Vector2d par(gsl_vector_get(x_par, 0), gsl_vector_get(x_par, 1));
  gsl_vector temp = grad_alpha_beta(par, w, quad_theta, post_w);
  gsl_vector_memcpy(f, &temp);
  return GSL_SUCCESS;
}

int Model::alpha_beta_df(const gsl_vector * x_par, void *params, gsl_matrix * J)
{
  const VectorXd w = ((struct rparams *) params) -> w;
  const VectorXd quad_theta = ((struct rparams *) params) -> quad_theta;
  const MatrixXd post_w = ((struct rparams *) params) -> post_w;
  const Vector2d par(gsl_vector_get(x_par, 0), gsl_vector_get(x_par, 1));
  gsl_matrix temp = jac_alpha_beta(par, w, quad_theta, post_w);
  gsl_matrix_memcpy(J, &temp);
  return GSL_SUCCESS;
}

int Model::alpha_beta_fdf(const gsl_vector * x_par, void *params, gsl_vector *
 f, gsl_matrix * J) {
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

gsl_vector Model::grad_a_b_gamma(VectorXd par, VectorXd x, VectorXd w, VectorXd
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

gsl_matrix Model::jac_a_b_gamma(VectorXd par, VectorXd x, VectorXd w, VectorXd
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

int Model::a_b_gamma_f(const gsl_vector * x_par, void *params, gsl_vector * f) {
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

int Model::a_b_gamma_df(const gsl_vector * x_par, void *params, gsl_matrix * J) 
{
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

int Model::a_b_gamma_fdf(const gsl_vector * x_par, void *params, gsl_vector * f,
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