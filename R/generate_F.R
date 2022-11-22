#'@export
generate_F_simple <- function(nblocks, type=c("nested", "difference",
                                              "checkers1", "checkers2")){

  type <- match.arg(type)
  if(type == "nested"){
    m <-  cbind(c(1, 1, 1), c(1, 1, 0), c(0, 0, 1))
  }else if(type == "difference"){
    m <- cbind(c(1, 1, 1), c(1, -1, 0), c(0, 0, 1))
  }else if(type == "checkers1"){
    m <- cbind(c(1, 1, 1),
               c(1, 1, 0 ),
               c(0, 1, 1))
  }else if(type == "checkers2"){
    m <-  cbind(c(1, 0, 1),
                c(1, 1, 0 ),
                c(0, 1, 1))
  }
  F_mat <- kronecker(diag(nblocks), m)
  return(F_mat)
}

#'@title Generate random F
#'@param g_F Function from which non-zero elements of F are generated
#'@param nz_factor Number of non-zero elements of each factor if F is to be generated.
#'@param omega Proportion of trait heritability explained by factors
#'@param h2_trait Trait heritability
#'@param pad Add single trait factors? (See details)
#'@details
#'Generate a random set of (at least) K factors for M traits.
#'The number of traits affected by each factor is given by \code{nz_factor}.
#'Each effect is chosen as a random draw from function g_F.
#'If any rows of the resulting matrix corresponding to non-zero elements of omega
#'are all zero and pad  = TRUE, single-trait factors are added.
#'Finally, the
#'matrix is re-scaled so that \code{colSums(F_mat^2) = omega*h2_trait}.
#'@return A matrix
#'@export
generate_random_F <- function(K, M, g_F= function(n){runif(n, -1, 1)},
                              nz_factor, omega, h2_trait,
                              pad = FALSE){
  nz_factor <- check_scalar_or_numeric(nz_factor, "nz_factor", K)
  omega <- check_scalar_or_numeric(omega, "omega", M)
  h2_trait <- check_scalar_or_numeric(h2_trait, "h2_trait", M)
  omega <- check_01(omega, "omega")
  h2_trait <- check_01(h2_trait, "h2_trait")
  srs <- omega*h2_trait
  F_mat <- sumstatFactors:::generate_F2(non_zero_by_factor = nz_factor,
                                        square_row_sums = srs,
                                        rfunc = g_F)
  if(any(rowSums(F_mat^2) == 0 & srs != 0) & pad){
    ix <- which(rowSums(F_mat^2)==0 & srs != 0)
    stF <- matrix(0, nrow = M, ncol = length(ix))
    for(i in seq_along(ix)) stF[ix[i],i] <- sign(srs[ix[i]])*sqrt(abs(srs[ix[i]]))
    F_mat <- cbind(F_mat, stF)
  }
  return(F_mat)
}

#'@export
generate_F2 <- function(non_zero_by_factor,
                        square_row_sums,
                        rfunc = function(n){runif(n, -1, 1)}){
  f <- function(n, nz){
    stopifnot(nz >= 1)
    x <- rep(0, n)
    ix <- sample(seq(n), size=nz, replace=F)
    x[ix] <- rfunc(nz)
    return(x)
  }
  ix <- which(square_row_sums > 0)
  #stopifnot(all(square_row_sums > 0))

  M <- length(square_row_sums)
  Mx <- length(ix)
  K <- length(non_zero_by_factor)

  F_mat <- matrix(0, M, K)
  F_matx <- sapply(seq(K), function(k){f(Mx, non_zero_by_factor[k])})
  F_mat[ix,] <- F_matx

  if(any(rowSums(F_mat !=0)==0)){
    missing_ix <- which(rowSums(F_mat !=0)==0)
    r2 <- rowSums(F_mat^2)
    F_mat <- F_mat*sqrt(square_row_sums)/sqrt(r2)
    F_mat[missing_ix,] <- 0
    return(F_mat)
  }
  r2 <- rowSums(F_mat^2)
  F_mat <- F_mat*sqrt(square_row_sums)/sqrt(r2)
  return(F_mat)
}
