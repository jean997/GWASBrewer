#'@title Generate beta hats from standardized or non-standardized direct SNP effects and LD
#'@param b_joint_std Matrix of standardized joint (causal) effects (dimension variants by traits)
#'@param b_joint  Matrix of non-standardized joint (causal) effects (dimension variants by traits). Supply only one of \code{b_joint} or
#'\code{b_joint_std}.
#'@param trait_corr Matrix of population trait correlation (traits by traits)
#'@param N Sample size, scalar, vector, or matrix. See \code{?sim_mv} for more details.
#'@param R_LD LD pattern (optional). See \code{?sim_mv} for more details.
#'@param af Allele frequencies (optional, allowed only if \code{R_LD} is missing). See \code{?sim_mv} for more details.
#'@details This function can be used to generate new GWAS results with the same effect sizes by passing in the \code{beta_joint} table
#'from a data set simulated using `sim_mv`. If the
#'original data are generated with af missing and no LD then the \code{beta_joint} table contains standardized effects. Otherwise
#'it contains non-standardized effects. Use the appropriate argument, either \code{b_joint_std} or \code{b_joint}.
#'@examples
#' # Use gen_bhat_from_b to generate new GWAS results with the same effect sizes.
#' N <- matrix(1000, nrow = 2, ncol =2)
#' G <- matrix(0, nrow = 2, ncol = 2)
#' R_E <- matrix(c(1, 0.8, 0.8, 1), nrow = 2, ncol = 2)
#' # original data
#' dat <- sim_mv(N = N, J = 20000, h2 = c(0.4, 0.3), pi = 1000/20000,
#'                G = G, R_E = R_E)
#' # data for second GWAS
#' # Since we didn't supply af or an LD pattern in the original GWAS,
#' # we have standardized effects.
#' dat_new <- gen_bhat_from_b(b_joint_std = dat$beta_joint, N = 40000,
#'                            trait_corr = dat$trait_corr)
#'@export
gen_bhat_from_b <- function(b_joint_std, b_joint,
                            trait_corr,  N,
                            R_LD = NULL,
                            af = NULL,
                            est_s = FALSE,
                            L_mat_joint_std = NULL,
                            theta_joint_std = NULL){

  #type <- match.arg(type, type)
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

  nn <- check_N(N, M)
  trait_corr <- check_matrix(trait_corr, "trait_corr", M, M)
  trait_corr <- check_psd(trait_corr, "trait_corr")


  # Correlation due to sample overlap is here
  R <- nn$Nc*trait_corr
  # sampling error on z-score without LD
  E_Z <- MASS::mvrnorm(n=J, mu = rep(0, M), Sigma = R)

  if(is.null(R_LD)){
    # compute sqrt(var(genos))
    af <- check_af(af, J)
    if(is.null(af)){
      sx <- rep(1, J)
    }else{
      sx <- sqrt(2*af*(1-af))
    }
    # se(beta_hat)_{ij} \approx 1/(sqrt(N_j)*sd(x_i))
    se_beta_hat <- kronecker(matrix(1/sx), matrix(1/sqrt(nn$N), nrow = 1))
    if(b_type == "non_std"){
      Z <- b_joint/se_beta_hat
    }else{
      Z <- t(t(b_joint_std)*sqrt(nn$N))
    }
    beta_hat <- (Z + E_Z)*se_beta_hat
    if(!is.null(af)){
      snp_info = data.frame(SNP = seq(J), AF = af)
    }else{
      snp_info = data.frame(SNP = seq(J), AF = NA)
    }
    ret <- list(beta_hat =beta_hat,
                se_beta_hat = se_beta_hat,
                sx = sx,
                R=R,
                Z = Z,
                snp_info = snp_info)
                #true_h2 = true_h2)

    if(est_s){
      ret$s_estimate <- estimate_s(N = N, beta_hat = beta_hat,
                                   trait_corr = trait_corr,
                                   af = af)
    }

    if(!is.null(L_mat_joint_std)){
      ret$L_mat <-(1/sx)*L_mat_joint_std
    }
    if(!is.null(theta_joint_std)){
      ret$theta <- (1/sx)*theta_joint_std
    }
    return(ret)
  }

  #If LD,

  if(is.null(af)) stop("Please provide af to go with R_LD.")

  ld_mat <- check_R_LD(R_LD, "matrix")
  ld_sqrt <- check_R_LD(R_LD, "sqrt")
  l <- check_R_LD(R_LD, "l")

  af <- check_af(af, sum(l), function_ok = FALSE)
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

  sx <- with(snp_info_full, sqrt(2*AF*(1-AF)))
  se_beta_hat <- kronecker(matrix(1/sx), matrix(1/sqrt(nn$N), nrow = 1)) # J by M


  # E_LD_Z are error terms which should have distribution
  # E_LD_Z ~ N(0, R)
  # Achieved by transforming E_Z_trait ~ N(0, 1) by
  # E_LD_Z = sqrt(R) E_Z with sqrt(R) = U D^{1/2}
  E_LD_Z <- lapply(seq(nb), function(i){
    tcrossprod(ld_sqrt[[block_info$block_index[i]]], t(E_Z[start_ix[i]:end_ix[i], ]))
  }) %>% do.call( rbind, .)


  if(b_type == "non_std"){
    Z <- b_joint/se_beta_hat
  }else{
    Z <- t(t(b_joint_std)*sqrt(nn$N))
  }

  # E[Z_marg] = R Z_joint
  Zm <- lapply(seq(nb), function(i){
          tcrossprod(ld_mat[[block_info$block_index[i]]], t(Z[start_ix[i]:end_ix[i], ]))
    }) %>% do.call( rbind, .)

  Z_hat <- Zm + E_LD_Z
  beta_hat <- Z_hat*se_beta_hat

  ret <- list(beta_hat =beta_hat,
              se_beta_hat = se_beta_hat,
              sx = sx,
              R=R,
              Z = Zm,
              snp_info = snp_info_full)

  if(est_s){
    ret$s_estimate <- estimate_s(N = N, beta_hat = beta_hat,
                                 trait_corr = trait_corr,
                                 R_LD = R_LD, af = af)
  }

  if(!is.null(L_mat_joint_std)){
    L_mat <- L_mat_joint_std*sx# S^-inv L (the N will cancel)
    L_mat <- lapply(seq(nb), function(i){
                     tcrossprod(ld_mat[[block_info$block_index[i]]], t(L_mat[start_ix[i]:end_ix[i], ]))
                      }) %>%
             do.call( rbind, .)
    L_mat <- L_mat/sx
    ret$L_mat <- L_mat
  }
  if(!is.null(theta_joint_std)){
    theta <- theta_joint_std*sx
    theta <- lapply(seq(nb), function(i){
                      tcrossprod(ld_mat[[block_info$block_index[i]]], t(theta[start_ix[i]:end_ix[i], ]))
                    }) %>%
      do.call( rbind, .)
    theta <- theta/sx
    ret$theta <- theta
  }

  return(ret)

}
