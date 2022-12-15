#'@export
gen_genos_mvn <- function(n, J, R_LD, af){

  af <- check_scalar_or_numeric(af, "af", J)
  af <- check_01(af, "af")

  if(is.null(R_LD)){
    X <- replicate(n = n, rbinom(n = J, size = 2, prob = af))
    return(list(X = X, af = af))
  }

  l <- purrr::map(R_LD, function(hh){nrow(hh$cov)}) %>% unlist()
  hdat <- R_LD_to_haplodat(R_LD, snp_info = data.frame(SNP = 1:J, AF = af))
  block_info <- assign_ld_blocks(l, J)
  nb <- length(l)
  if(!is.null(block_info$last_block_info)){
    b <- block_info$last_block_info[1]
    x <- block_info$last_block_info[2]
    last_block <- hdat[[b]]$cor[seq(x), seq(x)]
    last_info <- data.frame(SNP = paste0("last", seq(x)), AF = hdat[[b]]$freqs[seq(x)])
    hdat_last <- R_LD_to_haplodat(R_LD = list(last_block), last_info)
    hdat[[nb + 1]] <- hdat_last[[1]]
    l <- c(l, x)
    nb <- nb + 1
  }

  X <- purrr::map_dfr(block_info$block_index, function(bi){
    xx <- hapsim_simple(n = 2*n, hap = hdat[[bi]])
    a <- xx[seq(1, 2*n, by = 2), ] + xx[seq(2, 2*n, by = 2), ]
    a <- data.frame(t(a))
    names(a) <- paste0("S", 1:n)
    return(a)
  })
  return(list(X = X, af = af[block_info$index]))


}


#'@title Generate GWAS data from standardized or non-standardized direct SNP effects and LD
#'@param b_joint_std Matrix of standardized joint (causal) effects (dimension variants by traits)
#'@param b_joint  Matrix of non-standardized joint (causal) effects (dimension variants by traits). Supply only one of \code{b_joint} or
#'\code{b_joint_std}.
#'@param trait_corr Matrix of population trait correlation (traits by traits)
#'@param N Sample size, scalar, vector, or special sample size format data frame, see details.
#'@param R_LD LD pattern (optional). See \code{?sim_mv} for more details.
#'@param snp_info (optional, required if \code{R_LD} is supplied).
#'@param af Allele frequencies (optional, allowed only if \code{R_LD} is missing). See \code{?sim_mv} for more details.
#'@details This function can be used to generate individual level GWAS data by passing in the \code{beta_joint} table
#'from a data set simulated using `sim_mv`. If the
#'original data are generated with af missing and no LD then the \code{beta_joint} table contains standardized effects. Otherwise
#'it contains non-standardized effects. Use the appropriate argument, either \code{b_joint_std} or \code{b_joint}.
#'@export
gen_gwas_from_b <- function(b_joint_std, b_joint,
                            N, V_E, R_E = NULL,
                            R_LD = NULL, snp_info = NULL, af = NULL,
                            sim_func = gen_genos_mvn,
                            calc_sumstats = TRUE){


  if(!missing(b_joint)){
    if(!missing(b_joint_std)){
      stop("Please provide only one of b_joint or b_joint_std")
    }
    b_joint <- check_matrix(b_joint, "b_joint")
    M <- ncol(b_joint)
    J <- nrow(b_joint)
    b_type <- "non_std"
  }else if(!missing(b_joint_std)){
    b_joint_std <- check_matrix(b_joint_std, "b_joint_std")
    M <- ncol(b_joint_std)
    J <- nrow(b_joint_std)
    b_type <- "std"
  }

  message(paste0("SNP effects provided for ", J, " SNPs and ", M, " traits."))

  nn <- check_N(N, M, allow_mat = FALSE)
  if(is.null(nn$Ndf)){
    nn$Ndf <- make_Ndf_indep(nn$N)
  }
  V_E <- check_scalar_or_numeric(V_E, "V_E", M)
  V_E <- check_01(V_E, "V_E")
  if(is.null(R_E)){
    R_E <- diag(M, nrow = M)
  }
  R_E <- check_matrix(R_E, "R_E", M, M)
  R_E <- check_psd(R_E, "R_E")

  V <- diag(sqrt(V_E), nrow = M) %*% R_E %*% diag(sqrt(V_E), nrow = M)
  eV <- eigen(V)

  ntotal <- sum(nn$Ndf$N)
  if(is.null(R_LD)){
    if(is.null(af)){
      stop("Must supply one of af or R_LD.\n")
    }
    if(class(af) == "function"){
      myaf <- af(J)
    }
    af <- check_scalar_or_numeric(af, "af", J)
    af <- check_01(af)
    X <- sim_func(ntotal, J, NULL, af)
  }else{
    if(is.null(snp_info)){
      stop("Please provide snp_info with R_LD.\n")
    }
    X <- sim_func(ntotal, J, R_LD, sim_info$AF)
  }
  genos <- check_matrix(X$X, "genos", ntotal, J)
  af <- check_scalar_or_numeric(X$af, "af", J)

  if(type == "std"){
    # Compute non-standardized effects
    sx <- sqrt(2*af*(1-af))
    b_joint <- b_joint_std/sx
  }

  study_id <- rep(1:nrow(nn$Ndf), nn$Ndf$N)
  study_ix <- lapply(1:M, function(m){
    ix <- which(nn$Ndf[paste0("trait_", m)] == TRUE)
    which(study_id %in% ix)
  })

  y_E <- matrix(rnorm(n = M*ntotal), nrow = ntotal)
  y_E <- y_E %*% diag(sqrt(eV$values), M) %*% t(eV$vectors)

  Y <- purrr::map_dfc(1:M, function(m){
    b <- b_joint[,m]
    y_gen <- t(t(X[study_ix[[m]],])*b) %>% rowSums()
    y <- y_gen + y_E[study_ix[[m]], m]
    my_y <- data.frame(y = rep(NA, ntotal))
    names(my_y)[1] <- paste0("y_", m)
    my_y[study_ix[[m]],1] <- y
    return(my_y)
  })


  R <- list(X = X, Y= Y, af = af)
  if(!calc_sumstats) return(R)

  ## Calculate GWAS estimates fast
  R$sumstats <- fast_lm(X, Y, check = FALSE)

  return(R)

}

#'@export
fast_lm <- function(X, Y, check = TRUE){
  if(check){
    X <- check_matrix(X, "X")
    N <- nrow(X)
    J <- ncol(X)
    Y <- check_matrix(Y, "Y", N)
    M <- ncol(Y)
  }else{
    N <- nrow(X)
    J <- ncol(X)
    M <- ncol(Y)
  }
  sumstats <- purrr::map_dfc(1:M, function(m){
    x_ix <- which(!is.na(Y[,m]))
    yc <- Y[x_ix,m] - mean(Y[x_ix,m])
    xc <- scale(X[x_ix, ], center = TRUE, scale = FALSE)
    xty <- colSums(xc*yc)
    xtx <- colSums(xc^2)
    bhat <- xty/xtx
    yhat <- t(t(xc)*bhat)
    r2 <- colSums((yhat - yc)^2)
    s2 <- r2/(xtx*(length(yc)-2))
    df <- data.frame(bhat = bhat, s = sqrt(s2))
    names(df) <- paste0(c("bhat_", "s_"), m)
    return(df)
  })
  return(sumstats)
}
