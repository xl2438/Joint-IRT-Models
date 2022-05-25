data {
  int<lower = 0> N;
  int<lower = 0> J;
  int x[N, J];
  int w[N, J];
}
parameters {
  row_vector[2] theta[N];
  vector<lower = -5, upper = 5>[J] b;
  vector<lower = -5, upper = 5>[J] ga;
  vector<lower = -5, upper = 5>[J] beta;
  corr_matrix[2] rho;
}
model {
  vector[2] mu = [0, 0]';
  vector[2] tau = [1, 1]';
  rho ~ lkj_corr(2);
  target += multi_normal_lpdf(theta | mu, quad_form_diag(rho, tau));
  for (i in 1:N) {
    for (j in 1:J) {
      target += bernoulli_logit_lpmf(w[i, j] | theta[i, 2] - ga[j]);
      target += bernoulli_logit_lpmf(x[i, j] | theta[i, 1] - b[j] + w[i, j] * beta[j]);
    }
  }
}
