#include<Rcpp.h>

// [[Rcpp::export]]
double test(double out) {
	return(out + 3);
}