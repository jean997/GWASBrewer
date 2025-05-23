#'@title Generate beta hats from standardized or non-standardized direct SNP effects and LD
#'@param b_joint Matrix of joint effect sizes
#'@param N Sample size, scalar, vector, matrix, or data.frame. See \code{?sim_mv} for more details.
#'@param trait_corr Matrix of population trait correlation (traits by traits)
#'@param R_LD LD pattern (optional). See \code{?sim_mv} for more details.
#'@param R_LD_eigen Optionally, also supply list of eigen decompositions. The program will not check that this matches R_LD.
#'@param af Allele frequencies (optional, allowed only if \code{R_LD} is missing). See \code{?sim_mv} for more details.
#'@param est_s Estimate standard errors?
#'@param input_geno_scale Genotype scale of effects in  b_joint
#'@param input_pheno_sd Phenotype sd for effects in b_joint
#'@param output_geno_scale Output genotype scale
#'@param output_pheno_sd Output phenotype sd
#'@details This has been made an internal function. To resample summary statistics use \code{resample_sumstats}.
gen_bhat_from_b <- function(b_joint,
                            N,
                            trait_corr = NULL,
                            R_LD = NULL,
                            R_LD_eigen = NULL,
                            af = NULL,
                            est_s = FALSE,
                            input_geno_scale = c("allele", "sd"),
                            input_pheno_sd = 1,
                            output_geno_scale = c("allele", "sd"),
                            output_pheno_sd = 1){

  input_geno_scale <- match.arg(input_geno_scale)
  output_geno_scale <- match.arg(output_geno_scale)

  if((input_geno_scale == "allele" | output_geno_scale == "allele") & is.null(af)){
    stop("af is required if input or output geno scales are allele.")
  }
  M <- ncol(b_joint)
  J <- nrow(b_joint)
  message(paste0("SNP effects provided for ", J, " SNPs and ", M, " traits."))

  input_pheno_sd <- check_scalar_or_numeric(input_pheno_sd, "input_pheno_sd", M)
  output_pheno_sd <- check_scalar_or_numeric(output_pheno_sd, "output_pheno_sd", M)

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

    A <- get_convert_matrix(input_geno_scale = input_geno_scale,
                            output_geno_scale = output_geno_scale,
                            input_pheno_sd = input_pheno_sd,
                            output_pheno_sd = output_pheno_sd,
                            J = J,
                            af = af)

    beta_joint <- A*b_joint

    if(output_geno_scale=="allele"){
      snp_info = data.frame(SNP = seq(J), AF = af)
      sx <- sqrt(2*af*(1-af))
    }else{
      snp_info  <- data.frame(SNP = seq(J), AF = NA)
      sx <- rep(1, J)
    }

    # se(beta_hat)_{ij} \approx sd(y_j)/(sqrt(N_j)*sd(x_i))
    if(length(nnz$nonzero_ix) != 0){
      se_beta_hat[,nnz$nonzero_ix] <- kronecker(1/sx, output_pheno_sd[nnz$nonzero_ix]/sqrt(nnz$N)) |> matrix(nrow = J, byrow = T)
      Z[,nnz$nonzero_ix] <- beta_joint[,nnz$nonzero_ix]/se_beta_hat[,nnz$nonzero_ix]
      beta_hat[, nnz$nonzero_ix] <- (Z[,nnz$nonzero_ix] + E_Z[,nnz$nonzero_ix])*se_beta_hat[,nnz$nonzero_ix]
      if(est_s){
        my_af <- NULL
        if(output_geno_scale == "allele") my_af <- af
        s_estimate[,nnz$nonzero_ix] <- estimate_s(N = nnz,
                                     beta_hat = row_times(beta_hat[,nnz$nonzero_ix], 1/output_pheno_sd[nnz$nonzero_ix]), # convert beta_hat to pheno_sd = 1 scale for estimate_s
                                     trait_corr = trait_corr[nnz$nonzero_ix,nnz$nonzero_ix],
                                     af = my_af)
        s_estimate <- row_times(s_estimate, output_pheno_sd)
      }
    }
    ret <- list(beta_hat =beta_hat,
                se_beta_hat = se_beta_hat,
                R=R,
                sx = sx,
                beta_joint = beta_joint,
                beta_marg = beta_joint,
                snp_info = snp_info,
                geno_scale = output_geno_scale,
                pheno_sd = output_pheno_sd)


    if(est_s){
      ret$s_estimate <- s_estimate
    }
    return(ret)
  }

  #If LD,

  if(is.null(af)) stop("Please provide af to go with R_LD.")

  ld_mat <- check_R_LD(R_LD, "matrix")
  if(!is.null(R_LD_eigen)){
    ld_sqrt <- check_R_LD(R_LD_eigen, "sqrt")
  }else{
    ld_sqrt <- check_R_LD(R_LD, "sqrt")
  }
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
    elb <- fast_eigen(as.matrix(last_block))
    ld_sqrt[[nblock + 1]] <- with(elb, vectors * rep(sqrt(values), each = nrow(vectors)))
  }

  snp_info <- data.frame(SNP = 1:sum(l), AF = af)
  snp_info$block <- rep(seq(length(l)), l)
  snp_info$ix_in_block <- sapply(l, function(nl){seq(nl)}) %>% unlist()
  snp_info <- snp_info[block_info$index,]
  snp_info$rep <- rep(block_info$rep, block_info$l)
  snp_info$SNP <- with(snp_info, paste0(SNP, ".", rep))

  nb <- length(block_info$l)
  start_ix <- cumsum(c(1, block_info$l[-nb]))
  end_ix <- start_ix + block_info$l-1

  # convert b_joint to the sd/pheno_sd = 1 scale
  Astd <- get_convert_matrix(input_geno_scale = input_geno_scale,
                          output_geno_scale = "sd",
                          input_pheno_sd = input_pheno_sd,
                          output_pheno_sd = 1,
                          J = J,
                          af = snp_info$AF)
  b_joint <- Astd*b_joint

  # beta_marg_std = R beta_joint (S R S^{-1} = R since S is just diag(1/sqrt(N)))
  b_marg <- lapply(seq(nb), function(i){
    tcrossprod(ld_mat[[block_info$block_index[i]]], t(b_joint[start_ix[i]:end_ix[i], ,drop =F]))
  }) %>% do.call( rbind, .)

  # convert to final scale
  Aunstd <- get_convert_matrix(input_geno_scale = "sd",
                               output_geno_scale = output_geno_scale,
                               input_pheno_sd = 1,
                               output_pheno_sd = output_pheno_sd,
                               J = J,
                               af = snp_info$AF)

  beta_marg <- Aunstd*b_marg
  beta_joint <- Aunstd*b_joint
  if(output_geno_scale == "allele"){
    sx <- with(snp_info, sqrt(2*AF*(1-AF)))
  }else{
    sx <- rep(1, J)
  }
  E_LD_Z <- Zm <- matrix(nrow = J, ncol = M)
  if(length(nnz$nonzero_ix) > 0){
    se_beta_hat[,nnz$nonzero_ix] <- kronecker(1/sx, output_pheno_sd[nnz$nonzero_ix]/sqrt(nnz$N)) |> matrix(nrow = J, byrow = T)
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
      s_estimate[,nnz$nonzero_ix] <- estimate_s(N = nnz,
                                                beta_hat = row_times(beta_hat[,nnz$nonzero_ix], 1/output_pheno_sd[nnz$nonzero_ix]), # convert beta_hat to pheno_sd = 1 scale for estimate_s
                                                trait_corr = trait_corr[nnz$nonzero_ix,nnz$nonzero_ix],
                                                R_LD = R_LD,
                                                af = af)
      s_estimate <- row_times(s_estimate, output_pheno_sd)
      if(output_geno_scale == "sd"){
        s_estimate <- s_estimate*with(snp_info, sqrt(2*AF*(1-AF)))
      }
    }
  }
  ret <- list(beta_hat =beta_hat,
              se_beta_hat = se_beta_hat,
              sx = sx,
              R=R,
              beta_joint = beta_joint,
              beta_marg = beta_marg,
              snp_info = snp_info,
              geno_scale = output_geno_scale,
              pheno_sd = output_pheno_sd)

  if(est_s){
    ret$s_estimate <- s_estimate
  }
  return(ret)
}

compute_R_times_mat <- function(R_LD, af, J, X){
  ld_mat <- check_R_LD(R_LD, "matrix")
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
  }

  snp_info <- data.frame(SNP = 1:sum(l), AF = af)
  snp_info$block <- rep(seq(length(l)), l)
  snp_info$ix_in_block <- sapply(l, function(nl){seq(nl)}) %>% unlist()
  snp_info <- snp_info[block_info$index,]
  snp_info$rep <- rep(block_info$rep, block_info$l)
  snp_info$SNP <- with(snp_info, paste0(SNP, ".", rep))

  nb <- length(block_info$l)
  start_ix <- cumsum(c(1, block_info$l[-nb]))
  end_ix <- start_ix + block_info$l-1

  X_marg <- lapply(seq(nb), function(i){
    tcrossprod(ld_mat[[block_info$block_index[i]]], t(X[start_ix[i]:end_ix[i], ,drop =F]))
  }) %>% do.call( rbind, .)
  return(X_marg)
}
