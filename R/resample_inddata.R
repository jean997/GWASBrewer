
#'@title Sample individual level data with joint effects matching a sim_mv object
#'@param sim_dat An object of class \code{sim_mv} (produced by \code{sim_mv} or \code{gen_bhat_from_b}).
#'@param N Sample size, scalar, vector, or special sample size format data frame, see details.
#'@param V_E Vector with length equal to the number of traits giving the environmental variance of each trait.
#'@param R_E Environmental correlation matrix, (traits by traits). If missing, R_E is assumed to be the identity.
#'@param R_LD LD pattern (optional). See \code{?sim_mv} for more details.
#'@param af Allele frequencies. This can be a scalar, vector or a function. For this function, af must always be
#'supplied.
#'@details This function can be used to generate individual level GWAS data by passing in the \code{beta_joint} table
#'from a data set simulated using `sim_mv`. If the
#'original data are generated with af missing and no LD then the \code{beta_joint} table contains standardized effects. Otherwise
#'it contains non-standardized effects. Use the appropriate argument, either \code{b_joint_std} or \code{b_joint}.
#'@examples
#' # Use gen_gwas_from_b to generate individual level data with given effect size.
#' Ndf <- data.frame(trait_1 = 1, trait_2 = 1, N = 10000)
#' G <- matrix(0, nrow = 2, ncol = 2)
#' R_E <- matrix(c(1, 0.8, 0.8, 1), nrow = 2, ncol = 2)
#' # original data
#' dat <- sim_mv(N = Ndf, J = 2000, h2 = c(0.4, 0.3), pi = 100/2000,
#'                G = G, R_E = R_E, af = function(n){rbeta(n, 1, 5)})
#'# Now generate GWAS data
#'gw_dat <- resample_inddata(dat, N = Ndf, calc_sumstats = TRUE)
#'@export
resample_inddata <- function(sim_dat, N,
                            R_LD = NULL, af = NULL,
                            sim_func = gen_genos_mvn,
                            calc_sumstats = FALSE){

  if(any(is.na(sim_dat$snp_info))){
    # sim_dat contains standardized effects
    b_joint_std <- sim_dat$beta_joint
    M <- ncol(b_joint_std)
    J <- nrow(b_joint_std)
    b_type <- "std"
  }else{
    b_joint <- sim_dat$beta_joint
    M <- ncol(b_joint)
    J <- nrow(b_joint)
    b_type <- "non_std"
  }

  message(paste0("SNP effects provided for ", J, " SNPs and ", M, " traits."))

  nn <- check_N(N, M, allow_mat = FALSE)
  if(is.null(nn$Ndf)){
    nn$Ndf <- make_Ndf_indep(nn$N)
  }
  V <- sim_dat$Sigma_E
  eV <- eigen(V)

  ntotal <- sum(nn$Ndf$N)
  if(is.null(R_LD)){
    if(is.null(af)){
      stop("Must supply allele frequencies (af).\n")
    }
    af <- check_af(af, J)
    X <- sim_func(ntotal, J, NULL, af)
  }else{
    if(is.null(af)){
      stop("Please provide af with R_LD.\n")
    }
    l <- check_R_LD(R_LD, "l")
    af <- check_af(af, sum(l), function_ok = FALSE)
    X <- sim_func(ntotal, J, R_LD, af)
  }
  af <- check_af(X$af, J, function_ok = FALSE)
  X <- check_matrix(X$X, "genos", ntotal, J)


  if(b_type == "std"){
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
  sumstats <- fast_lm(X, Y, check = FALSE)
  R$beta_hat <- dplyr::select(sumstats, all_of(paste0("bhat_", 1:M)) ) %>% as.matrix()
  R$se_beta_hat <- dplyr::select(sumstats, all_of(paste0("s_", 1:M)) ) %>% as.matrix()
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


#'@export
gen_genos_mvn <- function(n, J, R_LD, af){

  if(is.null(R_LD)){
    af <- check_af(af, J)
    X <- replicate(n = n, rbinom(n = J, size = 2, prob = af)) %>% t()
    if(J == 1) X <- t(X)
    return(list(X = X, af = af))
  }
  # Checks happen in R_LD_to_haplodat
  l <- check_R_LD(R_LD, "l")
  #af <- check_af(af, sum(l), function_ok = FALSE)


  hdat <- R_LD_to_haplodat(R_LD, af = af)
  block_info <- assign_ld_blocks(l, J)
  nb <- length(l)
  if(!is.null(block_info$last_block_info)){
    b <- block_info$last_block_info[1]
    x <- block_info$last_block_info[2]
    last_block <- hdat[[b]]$cor[seq(x), seq(x)]
    last_af<- hdat[[b]]$freqs[seq(x)]
    hdat_last <- R_LD_to_haplodat(R_LD = list(last_block), last_af)
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

  return(list(X = t(X), af = af[block_info$index]))


}
