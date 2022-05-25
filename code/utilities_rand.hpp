#ifndef UTILITIES_RAND_HPP
#define UTILITIES_RAND_HPP
#include <random>

double rlnorm(double mu, double s) {
  std::default_random_engine generator;
  std::lognormal_distribution<double> distribution(mu, s);
  return distribution(generator);
}

double rnorm(double mu, double s) {
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(mu, s);
  return distribution(generator);
}


#endif