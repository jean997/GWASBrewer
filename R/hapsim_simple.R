
#' Generate haplotype data from a haplotype object
#' @description Generates haplotype data from a haplotype object.
#' @param n The number of haplotypes to generate.
#' @param hap A haplotype object containing allele frequencies and covariance matrix.
#' @param seed An optional seed for reproducibility.
#' @return A matrix of haplotype data, where rows correspond to haplotypes and columns to loci.
#' @details The function generates haplotype data based on the allele frequencies and covariance matrix provided in the `hap` object. It uses a multivariate normal distribution to generate the haplotypes.
#' @examples
#' # Example usage:
#' hap <- list(freqs = c(0.1, 0.2), cov = matrix(c(1, 0.5, 0.5, 1), nrow = 2))
#' hapsim_simple(n = 10, hap = hap, seed = 123)
hapsim_simple <- function(n, hap, seed = NULL){
  if(!is.null(seed))
    set.seed(seed)
  nloci <- length(hap$freqs)
  if(is.null(hap$eCov)){
    hap$eCov <- fast_eigen(hap$cov)
  }
  A <- matrix(stats::rnorm(n*nloci), nrow = n)
  sqrt_eV <- sqrt(pmax(hap$eCov$values, 0))
  A <- udvt(hap$eCov$vectors, sqrt_eV, A)
  A <- t(A)
  quants <- stats::qnorm(hap$freqs)
  Y <- t( t(A) > quants)
  return(Y)
}

#' Convert LD matrices to haplotype data
#' @description Converts a list of LD matrices to haplotype data.
#' @param R_LD A list of LD matrices, each matrix corresponding to a set of loci.
#' @param af A vector of allele frequencies corresponding to the loci in `R_LD`.  
#' @return A list of haplotype data, where each element corresponds to a set of loci.
#' @details The function formats LD matrices into hapsim formatted input data. 
R_LD_to_haplodat <- function(R_LD, af){
  ld_mat <- check_R_LD(R_LD, "matrix")
  l <- check_R_LD(R_LD, "l")
  af <- check_af(af, sum(l), function_ok = FALSE)

  start <- cumsum(c(1, l))[-(length(l) + 1)]
  stop <- start + l -1
  hdat <- lapply(seq(length(ld_mat)), function(i){
    nloci <- l[i]
    C <- ld_mat[[i]]
    P <- 1-af[start[i]:stop[i]]
    Q <- stats::qnorm(P)
    null.mat <- matrix(0, nrow = nloci, ncol = nloci)
    vmat <- .C("covariance", as.integer(nloci), as.double(C),
               as.double(P), as.double(Q), rlt = as.double(null.mat),
               PACKAGE = "hapsim")$rlt
    V <- matrix(vmat, nrow = nloci, ncol = nloci)
    if (!hapsim::checkpd(V)){
      warning("Coercing LD matrix to feasible values.\n")
      V <- hapsim::makepd(V)
    }
    eV <- fast_eigen(V)
    return(list(freqs  = P, cor = C, cov = V, eCov = eV))
  })
  return(hdat)
}
