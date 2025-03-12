#' Fast Eigen Decomposition using RcppArmadillo
#'
#' Computes the full eigen decomposition of a symmetric matrix using C++.
#' This function is faster than base R's `eigen()`.
#' This function will not check that your matrix is numeric or symmetric.
#'
#' @param m A symmetric numeric matrix.
#' @return A list with two components: `values` (eigenvalues) and `vectors` (eigenvectors).
#' @export
fast_eigen <- function(m) {
  #if (!is.matrix(m) || !is.numeric(m)) stop("Input must be a numeric matrix")
  #if (!all(abs(m - t(m)) < 1e-10)) stop("Matrix must be symmetric")
  x <- .Call(`_GWASBrewer_fast_eigen`, m)
  class(x) <- "eigen"
  return(x)
}

#' Fast Eigenvalues using RcppArmadillo
#'
#' Computes eigen values only of a symmetric matrix using C++.
#' This function will not check that your matrix is numeric or symmetric.
#'
#' @param m A symmetric numeric matrix.
#' @return A vector of eigen values
#' @export
fast_eigen_vals <- function(m) {
  .Call(`_GWASBrewer_fast_eigen_vals`, m)
}

# computes U %*% diag(d) %*% t(V)
udvt <- function(U, d, V){
  tcrossprod((U * rep(d, each = nrow(U))), V)
}

# computes diag(s) %*% R %*% diag(x)
sRx <- function(s, R, x){
  t( t(R*s) *x)
}