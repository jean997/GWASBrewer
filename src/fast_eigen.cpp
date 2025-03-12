#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List fast_eigen(const arma::mat& X) {
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, X);
  
  return Rcpp::List::create(Rcpp::Named("values") = eigval,
                            Rcpp::Named("vectors") = eigvec);
}
