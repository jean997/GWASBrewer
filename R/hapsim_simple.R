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
