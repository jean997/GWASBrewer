hapsim_simple <- function(n, hap, seed = NULL){
  if(!is.null(seed))
    set.seed(seed)
  nloci <- length(hap$freqs)
  if(is.null(hap$eCov)){
    hap$eCov <- eigen(hap$cov)
  }
  A <- matrix(rnorm(n*nloci), nrow = n)
  A <- hap$eCov$vectors %*% diag(sqrt(pmax(hap$eCov$values, 0)), nloci) %*% t(A)
  A <- t(A)
  quants <- qnorm(hap$freqs)
  Y <- t( t(A) > quants)
  return(Y)
}


R_LD_to_haplodat <- function(R_LD, snp_info){
  ld_mat <- check_R_LD(R_LD, "matrix")
  l <- sapply(ld_mat, nrow)
  snp_info <- check_snpinfo(snp_info, l)

  start <- cumsum(c(1, l))[-(length(l) + 1)]
  stop <- start + l -1
  hdat <- lapply(seq(length(ld_mat)), function(i){
    nloci <- l[i]
    C <- ld_mat[[i]]
    P <- 1-snp_info$AF[start[i]:stop[i]]
    Q <- qnorm(P)
    null.mat <- matrix(0, nrow = nloci, ncol = nloci)
    vmat <- .C("covariance", as.integer(nloci), as.double(C),
               as.double(P), as.double(Q), rlt = as.double(null.mat),
               PACKAGE = "hapsim")$rlt
    V <- matrix(vmat, nrow = nloci, ncol = nloci)
    if (!hapsim:::checkpd(V))
      V <- hapsim:::makepd(V)
    eV <- eigen(V)
    return(list(freqs  = P, cor = C, cov = V, eCov = eV))
  })
  return(hdat)
}
