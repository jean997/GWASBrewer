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
  stopifnot(all(square_row_sums > 0))

  M <- length(square_row_sums)
  K <- length(non_zero_by_factor)


  F_mat <- sapply(seq(K), function(k){f(M, non_zero_by_factor[k])})

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
