#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec fast_eigen_vals(const arma::mat& X) {
  return arma::eig_sym(X);  // Compute only eigenvalues
}