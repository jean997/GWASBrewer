#'@title Generate beta hats from standardized or non-standardized direct SNP effects and LD
#'@param b_joint_std Matrix of standardized joint (causal) effects (dimension variants by traits)
#'@param b_joint  Matrix of non-standardized joint (causal) effects (dimension variants by traits). Supply only one of \code{b_joint} or
#'\code{b_joint_std}.
#'@param trait_corr Matrix of population trait correlation (traits by traits)
#'@param N Sample size, scalar, vector, matrix, or data.frame. See \code{?sim_mv} for more details.
#'@param R_LD LD pattern (optional). See \code{?sim_mv} for more details.
#'@param af Allele frequencies (optional, allowed only if \code{R_LD} is missing). See \code{?sim_mv} for more details.
#'@details This function can be used to generate new GWAS results with the same effect sizes by passing in the \code{beta_joint} table
#'from a data set simulated using `sim_mv`. If the
#'original data are generated with af missing and no LD then the \code{beta_joint} table contains standardized effects. Otherwise
#'it contains non-standardized effects. Use the appropriate argument, either \code{b_joint_std} or \code{b_joint}.
#'
#'Future version: only accepts standardized effects. Returns effects in units of either per allele or per genotype sd. Phenotype unit is always 1.
#'@examples
#' # Use gen_bhat_from_b to generate new GWAS results with the same effect sizes.
#' N <- matrix(1000, nrow = 2, ncol =2)
#' R_E <- matrix(c(1, 0.8, 0.8, 1), nrow = 2, ncol = 2)
#' # original data
#' dat <- sim_mv(G = 2,
#'               N = N,
#'               J = 20000,
#'               h2 = c(0.4, 0.3),
#'               pi = 1000/20000,
#'               R_E = R_E)
#' # data for second GWAS
#' # Since we didn't supply af or an LD pattern in the original GWAS,
#' # we have standardized effects.
#' dat_new <- gen_bhat_from_b(b_joint_std = dat$beta_joint, N = 40000,
#'                            trait_corr = dat$trait_corr)
gen_bhat_from_b <- function(b_joint,
                            N,
                            trait_corr = NULL,
                            R_LD = NULL,
                            af = NULL,
                            est_s = FALSE,
                            L_mat_joint_std = NULL,
                            theta_joint_std = NULL,
                            return_geno_unit = c("allele", "sd"),
                            return_pheno_sd = 1){

  return_geno_unit <- match.arg(return_geno_unit)
  M <- ncol(b_joint)
  J <- nrow(b_joint)
  message(paste0("SNP effects provided for ", J, " SNPs and ", M, " traits."))

  nn <- check_N(N, M)
  nnz <- subset_N_nonzero(nn)

  ## Check for sample overlap
  if(!nn$overlap){
    trait_corr <- diag(M)
  }else{
    if(is.null(trait_corr)){
      stop("trait_corr required due to sample overlap.")
    }
  }
  trait_corr <- check_matrix(trait_corr, "trait_corr", M, M)
  trait_corr <- check_psd(trait_corr, "trait_corr")

  beta_hat <- se_beta_hat <- s_estimate <- Z <- E_Z <- matrix( nrow = J, ncol = M)
  R <- matrix(0, nrow = M, ncol = M)

  Mz <- length(nnz$nonzero_ix)
  # Correlation due to sample overlap
  R[nnz$nonzero_ix, nnz$nonzero_ix] <- Rz <- nnz$Nc*trait_corr[nnz$nonzero_ix, nnz$nonzero_ix]
  # sampling error on z-score without LD
  if(length(nnz$nonzero_ix) > 0) E_Z[,nnz$nonzero_ix] <- MASS::mvrnorm(n=J, mu = rep(0, Mz), Sigma = Rz)

  if(is.null(R_LD)){
    # compute sqrt(var(genos))
    af <- check_af(af, J)
    if(return_geno_unit == "sd"){ ## af is ignored if no LD and return unit is sd
      if(!is.null(af)){
        warning("af ignored if return_geno_unit = sd and there is no LD.\n")
      }
      sx <- rep(1, J)
      snp_info  <- data.frame(SNP = seq(J), AF = NA)
    }else{
      sx <- sqrt(2*af*(1-af))
      snp_info = data.frame(SNP = seq(J), AF = af)
    }
    beta_joint <- b_joint/sx
    # se(beta_hat)_{ij} \approx sd(y_j)/(sqrt(N_j)*sd(x_i))
    if(length(nnz$nonzero_ix) != 0){
      se_beta_hat[,nnz$nonzero_ix] <- kronecker(matrix(1/sx), matrix(1/sqrt(nnz$N), nrow = 1))
      Z[,nnz$nonzero_ix] <- beta_joint[,nnz$nonzero_ix]/se_beta_hat[,nnz$nonzero_ix]
      beta_hat[, nnz$nonzero_ix] <- (Z[,nnz$nonzero_ix] + E_Z[,nnz$nonzero_ix])*se_beta_hat[,nnz$nonzero_ix]
      if(est_s){
        my_af <- NULL
        if(return_geno_unit == "allele") my_af <- af
        s_estimate[,nnz$nonzero_ix] <- estimate_s(N = nnz,
                                     beta_hat = beta_hat[,nnz$nonzero_ix],
                                     trait_corr = trait_corr[nnz$nonzero_ix,nnz$nonzero_ix],
                                     af = my_af)

      }
    }
    ret <- list(beta_hat =beta_hat,
                se_beta_hat = se_beta_hat,
                R=R,
                sx = sx,
                beta_joint = beta_joint,
                beta_marg = beta_joint,
                snp_info = snp_info,
                geno_scale = return_geno_unit,
                pheno_sd = rep(1, M))


    if(est_s){
      ret$s_estimate <- s_estimate
    }
    if(!is.null(L_mat_joint_std)){
      if(return_geno_unit == "allele"){
        ret$L_mat <-(1/sx)*L_mat_joint_std
      }else{
        ret$L_mat <- L_mat_joint_std
      }
    }
    if(!is.null(theta_joint_std)){
      if(return_geno_unit == "allele"){
        ret$theta <- (1/sx)*theta_joint_std
      }else{
        ret$theta <- theta_joint_std
      }
    }
    return_pheno_sd <- check_scalar_or_numeric(return_pheno_sd, "return_pheno_sd", M)
    if(!all(return_pheno_sd == 1)){
      ret <- rescale_sumstats(ret,  output_geno_scale = return_geno_unit, output_pheno_sd = return_pheno_sd)
    }
    return(ret)
  }

  #If LD,

  if(is.null(af)) stop("Please provide af to go with R_LD.")

  ld_mat <- check_R_LD(R_LD, "matrix")
  ld_sqrt <- check_R_LD(R_LD, "sqrt")
  l <- check_R_LD(R_LD, "l")

  af <- check_af(af, sum(l), function_ok = FALSE)

  ## LD bookkeeping
  nblock <- length(l)
  block_info <- assign_ld_blocks(l, J)
  if(!is.null(block_info$last_block_info)){
    b <- block_info$last_block_info[1]
    x <- block_info$last_block_info[2]
    last_block <- ld_mat[[b]][seq(x), seq(x)]
    ld_mat[[nblock + 1]] <- last_block
    elb <- eigen(last_block)
    ld_sqrt[[nblock + 1]] <- with(elb, vectors %*% diag(sqrt(values)))
  }

  snp_info <- data.frame(SNP = 1:sum(l), AF = af)
  snp_info$block <- rep(seq(length(l)), l)
  snp_info$ix_in_block <- sapply(l, function(nl){seq(nl)}) %>% unlist()
  snp_info_full <- snp_info[block_info$index,]
  snp_info_full$rep <- rep(block_info$rep, block_info$l)
  snp_info_full$SNP <- with(snp_info_full, paste0(SNP, ".", rep))

  nb <- length(block_info$l)
  start_ix <- cumsum(c(1, block_info$l[-nb]))
  end_ix <- start_ix + block_info$l-1

  sx <- case_when(return_geno_unit== "allele" ~ with(snp_info_full, sqrt(2*AF*(1-AF))),
                  TRUE ~ rep(1, J))

  # beta_marg_std = R beta_joint (S R S^{-1} = R since S is just diag(1/sqrt(N)))
  b_marg <- lapply(seq(nb), function(i){
    tcrossprod(ld_mat[[block_info$block_index[i]]], t(b_joint[start_ix[i]:end_ix[i], ,drop =F]))
  }) %>% do.call( rbind, .)

  beta_marg <- b_marg/sx
  beta_joint <- b_joint/sx

  E_LD_Z <- Zm <- matrix(nrow = J, ncol = M)
  if(length(nnz$nonzero_ix) > 0){
    se_beta_hat[,nnz$nonzero_ix] <- kronecker(matrix(1/sx), matrix(1/sqrt(nnz$N), nrow = 1)) # J by M
    # E_LD_Z are error terms which should have distribution E_LD_Z ~ N(0, R)
    # Achieved by transforming E_Z_trait ~ N(0, 1) by
    # E_LD_Z = sqrt(R) E_Z with sqrt(R) = U D^{1/2}
    E_LD_Z[,nnz$nonzero_ix] <- lapply(seq(nb), function(i){
      tcrossprod(ld_sqrt[[block_info$block_index[i]]], t(E_Z[start_ix[i]:end_ix[i], nnz$nonzero_ix,drop = F ]))
    }) %>% do.call( rbind, .)

    Z[,nnz$nonzero_ix] <- beta_joint[,nnz$nonzero_ix]/se_beta_hat[,nnz$nonzero_ix]
    Zm[,nnz$nonzero_ix] <- beta_marg[,nnz$nonzero_ix]/se_beta_hat[,nnz$nonzero_ix]
    Z_hat <- Zm + E_LD_Z
    beta_hat[,nnz$nonzero_ix] <- Z_hat[,nnz$nonzero_ix]*se_beta_hat[,nnz$nonzero_ix]
    if(est_s){
      s_estimate[,nnz$nonzero_ix] <- estimate_s(N = nnz, beta_hat = beta_hat[,nnz$nonzero_ix],
                                   trait_corr = trait_corr[nnz$nonzero_ix,nnz$nonzero_ix],
                                   R_LD = R_LD, af = af)
      if(return_geno_unit == "sd"){
        s_estimate <- s_estimate*with(snp_info, sqrt(2*AF*(1-AF)))
      }
    }
  }
  ret <- list(beta_hat =beta_hat,
              se_beta_hat = se_beta_hat,
              sx = sx,
              R=R,
              beta_joint = b_joint/sx,
              beta_marg = b_marg/sx,
              snp_info = snp_info_full,
              geno_scale = return_geno_unit,
              pheno_sd = rep(1, M))

  if(est_s){
    ret$s_estimate <- s_estimate
  }
  if(!is.null(L_mat_joint_std)){
    L_mat <- L_mat_joint_std*sx# S^-inv L (the N will cancel)
    L_mat <- lapply(seq(nb), function(i){
          tcrossprod(ld_mat[[block_info$block_index[i]]], t(L_mat[start_ix[i]:end_ix[i], ,drop = F]))
    }) %>%
    do.call( rbind, .)
    L_mat <- L_mat/sx
    ret$L_mat <- L_mat
  }
  if(!is.null(theta_joint_std)){
    theta <- theta_joint_std*sx
    theta <- lapply(seq(nb), function(i){
        tcrossprod(ld_mat[[block_info$block_index[i]]], t(theta[start_ix[i]:end_ix[i], ,drop = F]))
      }) %>%
      do.call( rbind, .)
    theta <- theta/sx
    ret$theta <- theta
  }
  return_pheno_sd <- check_scalar_or_numeric(return_pheno_sd, "return_pheno_sd", M)
  if(!all(return_pheno_sd ==1)){
    ret <- rescale_sumstats(ret,  output_geno_scale = return_geno_unit, output_pheno_sd = return_pheno_sd)
  }
  return(ret)
}
