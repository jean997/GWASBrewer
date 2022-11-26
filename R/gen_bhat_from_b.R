#'@title Generate beta hats from standardized direct SNP effects and LD
#'@param b_joint_std variants by traits matrix of direct joint SNP effects.
#'@export
gen_bhat_from_b <- function(b_joint_std, b_joint,
                            trait_corr,  N,
                            R_LD = NULL, snp_info = NULL, af = NULL,
                            L_mat_joint_std = NULL, theta_joint_std = NULL){

  if(!missing(b_joint)){
    if(!missing(beta_joint_std)){
      stop("Please provide only one of b_joint or b_joint_std")
    }
    if(!"matrix" %in% class(b_joint)){
      stop("b_joint_std must have class matrix\n")
    }
    M <- ncol(b_joint)
    J <- nrow(b_joint)
    b_type <- "non_std"
  }else if(!missing(b_joint_std)){
    if(!"matrix" %in% class(b_joint_std)){
      stop("b_joint_std must have class matrix\n")
    }
    # Expected z-score
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
    if(is.null(af)){
      sx <- rep(1, J)
    }else{
      af <- check_scalar_or_numeric(af, "af", J)
      sx <- sqrt(2*af*(1-af))
      if(length(sx) == 1) sx <- rep(sx, J)
    }
    # se(beta_hat)_{ij} = 1/(sqrt(N_j)*sd(x_i))
    se_beta_hat <- kronecker(matrix(1/sx), matrix(1/sqrt(nn$N), nrow = 1))
    if(b_type == "non_std"){
      Z <- b_joint/se_beta_hat
    }else{
      Z <- t(t(b_joint_std)*sqrt(nn$N))
    }
    beta_hat <- (Z + E_Z)*se_beta_hat

    ret <- list(beta_hat =beta_hat,
                se_beta_hat = se_beta_hat,
                sx = sx,
                R=R,
                Z = Z)

    if(!is.null(L_mat_joint_std)){
      ret$L_mat <-(1/sx)*L_mat_joint_std
    }
    if(!is.null(theta_joint_std)){
      ret$theta <- (1/sx)*theta_joint_std
    }
    return(ret)
  }

  #If LD,

  if(is.null(snp_info)) stop("Please provide snp_info to go with R_LD.")
  R_LD <- check_R_LD(R_LD)
  l <- sapply(R_LD, function(e){length(e$values)})
  snp_info <- check_snpinfo(snp_info, l)
  nblock <- length(l)

  #Bookkeeping: Figure out how much/how many replicates of supplied LD we need
  ld_size <- sum(l)
  full_reps <- floor(J/ld_size) # Recall l is list of block sizes
  remainder <- J  - full_reps*ld_size

  block_index <- rep(seq(nblock), full_reps)
  if(remainder > 0){
    if(l[1] >= remainder){
      blocks_rem <- 0
      final_remainder <- remainder
    }else{
      blocks_rem <- max(which(cumsum(l) <= remainder)) # full blocks in last partial repeat
      final_remainder <- remainder-cumsum(l)[blocks_rem] # partial block in last partial repeat
      block_index <- c(block_index, seq(blocks_rem))
    }
    if(final_remainder > 0){
      block_index <- c(block_index, nblock + 1)
      last_block <- with(R_LD[[blocks_rem + 1]], (vectors %*% diag(values) %*% t(vectors))[1:final_remainder, 1:final_remainder])
      R_LD[[nblock + 1]] <- eigen(last_block)
      l <- c(l, final_remainder)[block_index] # l is now lengths of blocks in data
    }else{
      l <- l[block_index] # l is now lengths of blocks in data
    }
    #Create SNP info for full data set
    snp_info_full <- snp_info[c(rep(seq(ld_size), full_reps), seq(remainder)),]
    if(full_reps > 0) snp_info_full$rep <- c(rep(seq(full_reps), each = ld_size), rep(full_reps + 1, remainder))
  }else{
    block_index <- rep(seq(nblock), full_reps)
    l <- l[block_index] # l is now lengths of blocks in data
    #Create SNP info for full data set
    snp_info_full <- snp_info[rep(seq(ld_size), full_reps),]
    if(full_reps > 0) snp_info_full$rep <- rep(seq(full_reps), each = ld_size)
  }
  if(full_reps == 0){
    snp_info_full$rep <- 1
  }

  snp_info_full$SNP <- with(snp_info_full, paste0(SNP, ".", rep))
  sx <- with(snp_info_full, sqrt(2*AF*(1-AF)))
  se_beta_hat <- kronecker(matrix(1/sx), matrix(1/sqrt(nn$N), nrow = 1)) # J by M

  # Multiply errors by square root of LD matrix
  E_LD_Z <- lapply(seq_along(block_index), function(i){
    #with(R_LD[[block_index[i]]], vectors %*% diag(sqrt(values)) %*% E_Z[start_ix[i]:end_ix[i], ])
    with(R_LD[[block_index[i]]], tcrossprod(vectors, tcrossprod(t(E_Z[start_ix[i]:end_ix[i], ]), diag(sqrt(values)))))
  }) %>% do.call( rbind, .)

  if(b_type == "non_std"){
    Z <- b_joint/se_beta_hat
  }else{
    Z <- t(t(b_joint_std)*sqrt(nn$N))
  }

  # Transform Z by LD matrix
  Z <- lapply(seq_along(block_index), function(i){
    #with(R_LD[[block_index[i]]], vectors %*% diag(values) %*% t(vectors) %*% Z[start_ix[i]:end_ix[i], ])
    with(R_LD[[block_index[i]]], tcrossprod(tcrossprod(vectors, tcrossprod(vectors, diag(values))), t(Z[start_ix[i]:end_ix[i], ])))
  }) %>% do.call( rbind, .)

  Z_hat <- Z + E_LD_Z
  beta_hat <- Z_hat*se_beta_hat

  ret <- list(beta_hat =beta_hat,
              se_beta_hat = se_beta_hat,
              sx = sx,
              R=R,
              Z = Z,
              snp_info = snp_info_full)

  if(!is.null(L_mat_joint_std)){
    L_mat <- L_mat_joint_std*sx# S^-inv L (the N will cancel)
    L_mat <- lapply(seq_along(block_index),
                    function(i){with(R_LD[[block_index[i]]],
                                     #vectors %*% diag(values) %*% t(vectors) %*% L_mat[start_ix[i]:end_ix[i], ])
                                     tcrossprod(tcrossprod(vectors, tcrossprod(vectors, diag(values))), t(L_mat[start_ix[i]:end_ix[i], ])))
                              }) %>%
             do.call( rbind, .)
    L_mat <- L_mat/sx
    ret$L_mat <- L_mat
  }
  if(!is.null(theta_joint_std)){
    theta <- theta_joint_std*sx
    theta <- lapply(seq_along(block_index),
                    function(i){with(R_LD[[block_index[i]]],
                                     #vectors %*% diag(values) %*% t(vectors) %*% theta[start_ix[i]:end_ix[i], ])
                                     tcrossprod(tcrossprod(vectors, tcrossprod(vectors, diag(values))), t(theta[start_ix[i]:end_ix[i], ])))
                    }) %>%
      do.call( rbind, .)
    theta <- theta/sx
    ret$theta <- theta
  }

  return(ret)

}
