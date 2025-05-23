## scaling: requires that beta_hat is allele, pheno_sd = 1 if af supplied
## otherwise beta_hat should be sd, pheno_sd = 1
estimate_s <- function(N, beta_hat,
                       trait_corr,
                       R_LD = NULL,
                       af = NULL){

  beta_hat <- check_matrix(beta_hat)
  J <- nrow(beta_hat)
  M <- ncol(beta_hat)
  nn <- check_N(N, M)
  trait_corr <- check_matrix(trait_corr, "trait_corr", M, M)
  trait_corr <- check_psd(trait_corr, "trait_corr")

  # s2_num_s <- 2/nn$N
  # s2_num_cor <- nn$Nc*trait_corr^2
  s2_S <- diag(sqrt(2/nn$N), nrow = M) %*% (nn$Nc*(trait_corr^2)) %*% diag(sqrt(2/nn$N), nrow = M)
  s2_num <- MASS::mvrnorm(n = 1, mu = rep(1, M), Sigma = s2_S)
  s2_num <- s2_num*(nn$N/(nn$N-2))
  #cat(length(s2_num), "\n")
  #
  s2_denom <- MASS::mvrnorm(n = J, mu = rep(0, M), Sigma = nn$Nc)
  if(is.null(R_LD)){
    af <- check_af(af, J)
    if(is.null(af)){
      v2x <- rep(1, J)
      s2x <- rep(1, J)
    }else{
      v2x <- 2*af*(1-af)
      s2x <- v2x*(1-v2x) # Variance of S2x
    }
    s2_denom_s <- kronecker(matrix(sqrt(s2x), nrow = J), matrix(sqrt(nn$N), nrow = 1))

    s2_denom_final <-  kronecker(matrix(v2x, nrow = J), matrix(nn$N, nrow = 1)) + s2_denom*s2_denom_s

    s2_estimate <- t( t(1/s2_denom_final)*s2_num) - t(t(beta_hat^2)/(nn$N-2))
    iter <- 0
    while(any(s2_estimate < 0)){
      iter <- iter + 1
      cat("updating to fix negaive s2: ", iter, "\n")
      i <- which(apply(s2_estimate, 1, min) < 0)
      s2_denom[i,] <- MASS::mvrnorm(n = length(i), mu = rep(0, M), Sigma = nn$Nc)
      s2_denom_final <-  kronecker(matrix(v2x, nrow = J), matrix(nn$N, nrow = 1)) + s2_denom*s2_denom_s
      s2_estimate <- t( t(1/s2_denom_final)*s2_num) - t(t(beta_hat^2)/(nn$N-2))
    }
    return(sqrt(s2_estimate))
  }

  ld_mat <- check_R_LD(R_LD, "matrix")
  l <- check_R_LD(R_LD, "l")
  nblock <- length(l)
  af <- check_af(af, sum(l), function_ok = FALSE)


  start_ix <- cumsum(c(1, l[-nblock]))
  end_ix <- start_ix + l -1
  v <- sqrt(2*af*(1-af))
  s2x <- (v^2)*(1-(v^2))
  one_minus_mu <- 1 - 2*af

  #s_denom <- lapply(seq(nb), function(i){
  sx_Rsqrt <- lapply(seq(nblock), function(i){
    s <- start_ix[i]
    p <- end_ix[i]
    vvt <- kronecker(matrix(v[s:p], ncol=1), matrix(v[s:p], nrow =1))
    oot <- kronecker(matrix(one_minus_mu[s:p], ncol=1), matrix(one_minus_mu[s:p], nrow =1))

    sx_R <- cov2cor(ld_mat[[i]]*vvt*oot + (ld_mat[[i]]^2)*(vvt^2)) 
    esx_R <- fast_eigen(sx_R)
    with(esx_R,  vectors*rep(sqrt(values), each = nrow(vectors)))

  })


  block_info <- assign_ld_blocks(l, J)
  if(!is.null(block_info$last_block_info)){
    b <- block_info$last_block_info[1]
    x <- block_info$last_block_info[2]

    s <- start_ix[b]
    p <- s + x -1
    vvt <- kronecker(matrix(v[s:p], ncol=1), matrix(v[s:p], nrow =1))
    oot <- kronecker(matrix(one_minus_mu[s:p], ncol=1), matrix(one_minus_mu[s:p], nrow =1))
    sx_Rl <- stats::cov2cor(ld_mat[[b]][seq(x), seq(x)]*vvt*oot + (ld_mat[[b]][seq(x), seq(x)]^2)*(vvt^2)) 
    esx_Rl <- fast_eigen(sx_Rl)
    sx_Rsqrtl <- with(esx_Rl,  vectors*rep(sqrt(values), each = nrow(vectors)))
            

    sx_Rsqrt[[nblock + 1]] <- sx_Rsqrtl
  }

  AF <- af[block_info$index]

  nb <- length(block_info$l)
  start_ix <- cumsum(c(1, block_info$l[-nb]))
  end_ix <- start_ix + block_info$l-1
  
  s2_estimate <- lapply(seq(nb), function(i){
    
    d2 <- tcrossprod(sx_Rsqrt[[block_info$block_index[i]]], t(s2_denom[start_ix[i]:end_ix[i], ]))
    s2x_i <- s2x[block_info$index[start_ix[i]:end_ix[i]]] # Variance of S2x
    v2x_i <- (v^2)[block_info$index[start_ix[i]:end_ix[i]]]
    s2_denom_si <- kronecker(matrix(sqrt(s2x_i), nrow = block_info$l[i]), matrix(sqrt(nn$N), nrow = 1))
    
    s2_denom_finali <-  kronecker(matrix(v2x_i, nrow = block_info$l[i]), matrix(nn$N, nrow = 1)) + d2*s2_denom_si
    
    x <- t( t(1/s2_denom_finali)*s2_num) - t(t(beta_hat[start_ix[i]:end_ix[i],]^2)/(nn$N-2))
    iter <- 0
    while(any(x < 0)){
      iter <- iter + 1
      cat("updating to fix negaive s2: ", i, " ", iter, "\n")
      s2_denom[start_ix[i]:end_ix[i], ] <- MASS::mvrnorm(n = block_info$l[i], mu = rep(0, M), Sigma = nn$Nc)
      
      d2 <- tcrossprod(sx_Rsqrt[[block_info$block_index[i]]], t(s2_denom[start_ix[i]:end_ix[i], ]))
      s2_denom_finali <-  kronecker(matrix(v2x_i, nrow = block_info$l[i]), matrix(nn$N, nrow = 1)) + d2*s2_denom_si
      
      x <- t( t(1/s2_denom_finali)*s2_num) - t(t(beta_hat[start_ix[i]:end_ix[i],]^2)/(nn$N-2))
    }
    return(x)
  }) %>% do.call(rbind, .)
  
  # s2_denom2 <- lapply(seq(nb), function(i){
  #   tcrossprod(sx_Rsqrt[[block_info$block_index[i]]], t(s2_denom[start_ix[i]:end_ix[i], ]))
  # }) %>% do.call( rbind, .)
  # 
  # 
  # s2x <- s2x[block_info$index] # Variance of S2x
  # v2x <- (v^2)[block_info$index]
  # s2_denom_s <- kronecker(matrix(sqrt(s2x), nrow = J), matrix(sqrt(nn$N), nrow = 1))
  # 
  # s2_denom_final <-  kronecker(matrix(v2x, nrow = J), matrix(nn$N, nrow = 1)) + s2_denom2*s2_denom_s
  # 
  # s2_estimate <- t( t(1/s2_denom_final)*s2_num) - t(t(beta_hat^2)/(nn$N-2))
 
  return(sqrt(s2_estimate))
}
