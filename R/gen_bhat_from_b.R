#'@title Generate beta hats from standardized direct SNP effects and LD
#'@param b_joint_std variants by traits matrix of direct joint SNP effects.
#'@param R traits by traits correlation matrix. This gives the correlation of z-scores due to sample overlap and environmental correlation.
#'This is output as the R object from sim_mv.
#'@export
gen_bhat_from_b <- function(b_joint_std, R,  N,
                            R_LD = NULL, snp_info = NULL, maf = NULL,
                            L_mat_joint_std = NULL, theta_joint_std = NULL){

  if(!"matrix" %in% class(b_joint_std)){
    stop("b_joint_std must have class matrix\n")
  }
  if(!"matrix" %in% class(R)){
    stop("R must have class matrix\n")
  }

  M <- ncol(b_joint_std)
  J <- nrow(b_joint_std)
  message(paste0("SNP effects provided for ", J, " SNPs and ", M, " traits."))

  if(length(N) == 1) N <- rep(N, M)
  else if(length(N) != M){
    stop(paste0("Inconsistent arguments: N has length ", length(N), " but b_joint_std has ", M, " columns."))
  }
  if(nrow(R) != M | ncol(R) != M | !Matrix::isSymmetric(R) | !all(diag(R) == 1)){
    stop(paste0("R should be an ", M , " times ", M, " correlation matrix."))
  }


  E_Z <- MASS::mvrnorm(n=J, mu = rep(0, M), Sigma = R)
  Z <- t(t(b_joint_std)*sqrt(N))

  if(missing(R_LD) | is.null(R_LD)){

    # generate MAF if required and compute sqrt(var(genos))
    if(is.null(maf)){
      sx <- rep(1, J)
    }else if(class(maf) == "numeric"){
      stopifnot(length(maf) %in% c(1, J))
      sx <- sqrt(2*maf*(1-maf))
      if(length(sx) == 1) sx <- rep(sx, J)
    }


    se_beta_hat <- kronecker(matrix(1/sx), matrix(1/sqrt(N), nrow = 1))
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

  if(missing(snp_info) | is.null(snp_info)) stop("Please prvide snp_info to go with R_LD.")
  l <- sapply(R_LD, function(e){length(e$values)})
  nblock <- length(R_LD)

  stopifnot(nrow(snp_info) == sum(l))
  stopifnot(all(c("SNP", "AF") %in% names(snp_info)))
  snp_info$block <- rep(seq(length(l)), l)
  snp_info$ix_in_block <- sapply(l, function(ll){seq(ll)}) %>% unlist()

  #If LD, introduce correlation

  #Bookkeeping: Figure out how much/how many replicates of supplied LD we need
  ld_size <- sum(l)
  full_reps <- floor(J/ld_size) # Recall l is list of block sizes
  remainder <- J  - full_reps*ld_size
  blocks_rem <- max(which(cumsum(l) <= remainder)) # full blocks in last partial repeat
  final_remainder <- remainder-cumsum(l)[blocks_rem] # partial block in last partial repeat

  last_block <- with(R_LD[[blocks_rem + 1]], (vectors %*% diag(values) %*% t(vectors))[1:final_remainder, 1:final_remainder])
  R_LD[[nblock + 1]] <- eigen(last_block)
  block_index <- c(rep(seq(nblock), full_reps), seq(blocks_rem), nblock + 1)
  l <- c(l, final_remainder)[block_index] # l is now lengths of blocks in data
  start_ix <- cumsum(c(1, l[-length(l)]))
  end_ix <- start_ix + l-1

  #Create SNP info for full data set
  snp_info_full <- snp_info[c(rep(seq(ld_size), full_reps), seq(remainder)),]
  if(full_reps == 0){
    snp_info_full$rep <- rep(1, remainder)
  }else{
    snp_info_full$rep <- c(rep(seq(full_reps), each = ld_size), rep(full_reps + 1, remainder))
  }
  snp_info_full$SNP <- with(snp_info_full, paste0(SNP, ".", rep))
  sx <- with(snp_info_full, sqrt(2*AF*(1-AF)))


  # Multiply errors by square root of LD matrix
  E_LD_Z <- lapply(seq_along(block_index), function(i){
    with(R_LD[[block_index[i]]], vectors %*% sqrt(diag(values)) %*% E_Z[start_ix[i]:end_ix[i], ])
  }) %>% do.call( rbind, .)

  # Transform Z by LD matrix
  Z <- lapply(seq_along(block_index), function(i){
    with(R_LD[[block_index[i]]], vectors %*% diag(values) %*% t(vectors) %*% Z[start_ix[i]:end_ix[i], ])
  }) %>% do.call( rbind, .)

  Z_hat <- Z + E_LD_Z

  se_beta_hat <- kronecker(matrix(1/sx), matrix(1/sqrt(N), nrow = 1)) # J by M
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
                                     vectors %*% diag(values) %*% t(vectors) %*% L_mat[start_ix[i]:end_ix[i], ])
                              }) %>%
             do.call( rbind, .)
    L_mat <- L_mat/sx
    ret$L_mat <- L_mat
  }
  if(!is.null(theta_joint_std)){
    theta <- theta_joint_std*sx
    theta <- lapply(seq_along(block_index),
                    function(i){with(R_LD[[block_index[i]]],
                                     vectors %*% diag(values) %*% t(vectors) %*% theta[start_ix[i]:end_ix[i], ])
                    }) %>%
      do.call( rbind, .)
    theta <- theta/sx
    ret$theta <- theta
  }

  return(ret)

}
