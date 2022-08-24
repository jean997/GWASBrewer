
#See simulating_data.Rmd for explanation
#'@title Simulate summary stats
#'@description Simulate summary statistcs for fully overlapping GWAS with no LD
#'@param F_mat factor matrix M by K
#'@param N GWAS sample size
#'@param J Total number of SNPs to generate
#'@param h_2_trait Heritability of each trait. Length M vector.
#'@param omega Proportion of SNP heritability mediated by factors. Length M vector.
#'@param h_2_factor Heritability of each factor. Length K vector.
#'@param pi_L Proportion of non-zero elements in L_k. Length K factor
#'@param pi_theta Proportion of non-zero elements in theta. Scalar or length M vector.
#'@param R_E Environmental trait correlation not mediated by factors. M by M pd matrix.
#'@param maf can either be a scalar in which case the same maf is used for all SNPS, NA in which case SNPs
#'are assumed scaled to variance 1, a function that takes a number n and returns n values between 0 and 1, or a
#'vector of length equal to the number of SNPs. maf is ignored if LD is provided.
#'@param R_LD List of eigen decompositions of LD correlation matrices, may be missing.
#'@param snp_info If R_LD is provided, provide a data frame with columns "SNP" and "AF"
#'@param g_F Function from which non-zero elements of F are generated
#'@param nz_factor Number of non-zero elements of each factor if F is to be generated.
#'@param sporadic_pleiotropy Allow a single SNP to affect multiple factors (default TRUE).
#'@details
#'
#'If F_mat is not provided, it will be generated using the `generate_F2` function.
#'In this case g_F and nz_factor must be provided. All of the elements
#'in F are generated iid from a mixture of a point mass at 0 and g_F. The matrix is then
#'rescaled according to the constraints.
#'
#'
#'With this setup it is possible to specify a setting that is impossible. Usually this occurs if the heritability
#'of the factors is low but the heritability of the traits is high leading to a contradiction. Right now the function
#'will just return an error if that happens.
#'
#'
#'@export
sim_sumstats_lf <- function(F_mat, N, J, h_2_trait, omega, h_2_factor,
                            pi_L, pi_theta,
                            R_E, maf = NA, R_LD = NULL, snp_info = NULL,
                            g_F, nz_factor,
                            overlap_prop =1,
                            sporadic_pleiotropy = TRUE){

  #Check parameters
  if(!missing(F_mat)){
    stopifnot("matrix" %in% class(F_mat))
    M <- nrow(F_mat)
    K <- ncol(F_mat)
    cat("Factor structure provided for ", M, " traits and ", K, " factors.\n")
  }else{
    if(missing(g_F) | missing(nz_factor) ){
      stop("If F_mat is missing please supply g_F and nz_factor")
    }
    K <- length(nz_factor)
    M <- length(h_2_trait)
    cat("Will generate L and F with ", J, " markers, ", M, " traits, and ", K, " factors.\n")
  }
  stopifnot(length(h_2_trait) == M)
  stopifnot(all(h_2_trait <= 1 & h_2_trait >= 0))


  stopifnot(length(omega) == M)
  stopifnot(all(omega >= 0 & omega <= 1))
  stopifnot(length(pi_L) == K)
  stopifnot(all(pi_L <= 1 & pi_L > 0))

  if(length(pi_theta) == 1){
    pi_theta <- rep(pi_theta, M)
  }else{
    stopifnot(length(pi_theta)== M)
  }
  stopifnot(all(pi_theta >=0 & pi_theta <=1))
  if(any(omega < 1 & pi_theta == 0)){
    stop("Omega is less than 1 for trait ", i, " so pi_theta must be greater than 0 for that trait.")
  }

  #R_E
  if(overlap_prop > 0){
    if(missing(R_E) | is.null(R_E)){
      message("R_E not provided but overlap_prop > 0. Using R_E = diag(ntrait) for no environmental covariance.")
      R_E <- diag(M)
    }
    if(missing(h_2_factor))('h_2_factor must be provided if overlap_prop > 0.')
    stopifnot(nrow(R_E) == M & ncol(R_E) == M)
    stopifnot(Matrix::isSymmetric(R_E))
    R_E_eig <- eigen(R_E)
    stopifnot(all(R_E_eig$values >= 0))

    stopifnot(length(h_2_factor) == K)
    stopifnot(all(h_2_factor <= 1 & h_2_factor >= 0))
  }

  #N
  if(length(N) == 1) N <- rep(N, M)
    else stopifnot(length(N) == M)

  #maf
  if(missing(R_LD) | is.null(R_LD)){
    if(is.na(maf)){
      sx <- rep(1, J)
    }else if(class(maf) == "numeric"){
      stopifnot(length(maf) %in% c(1, J))
      sx <- sqrt(2*maf*(1-maf))
      if(length(sx) == 1) sx <- rep(sx, J)
    }else if(class(maf) == "function"){
      af <- maf(J)
      sx <- sqrt(2*af*(1-af))
    }
  }else{
    if(missing(snp_info) | is.null(snp_info)) stop("Please prvide snp_info to go with R_LD.")
    l <- sapply(R_LD, function(e){length(e$values)})
    stopifnot(nrow(snp_info) == sum(l))
    stopifnot(all(c("SNP", "AF") %in% names(snp_info)))
    snp_info$block <- rep(seq(length(l)), l)
  }


  #Re-scale F or generate it if it is missing
  if(missing(F_mat)){
    F_mat <- sumstatFactors:::generate_F2(non_zero_by_factor = nz_factor,
                         square_row_sums = omega*h_2_trait,
                         rfunc = g_F)
    if(any(rowSums(F_mat^2) == 0)){
      ix <- which(rowSums(F_mat^2)==0)
      omega[ix] <- 0
    }
  }else{
    if(any(rowSums(F_mat == 0) == K & omega >0)){
      stop("One row of F is zero but corresponds to non-zero omega\n")
    }
    #Re scale rows of F
    scale <- sqrt(omega*h_2_trait/rowSums(F_mat^2))
    scale[omega == 0] <- 0
    F_mat <- F_mat*scale
  }


  #Generate theta
  sigma_theta <- sqrt( (1/(pi_theta*J)) * (1-omega)*h_2_trait)
  sigma_theta[omega == 1] <- 0

  if(all(sigma_theta == 0)){
    theta <- matrix(0, nrow = J, ncol = M)
  }else if(sporadic_pleiotropy){
    theta <- purrr::map(seq(M), function(i){
      t <- rbinom(n=J, size=1, prob = pi_theta[i])
      n <- sum(t==1)
      t[t==1] <- rnorm(n=n, mean=0, sd = sigma_theta[i])
      return(t)
    }) %>% do.call(cbind, .)
  }else{
    ix <- which(sigma_theta > 0)
    if(sum(pi_theta[ix]) > 1){
      stop("You have requested too many traits and too many causal variants to use sporadic_pleiotrop = FALSE.\n")
    }
    p <- sum(pi_theta[ix])
    nz_theta_ix <- sample(c(0, ix), size = J, replace = TRUE, prob = c(1-p, pi_theta[ix]) )
    theta <- purrr::map(seq(M), function(i){
      t <- which(nz_theta_ix == i)
      n <- length(t)
      val <- rep(0, J)
      if(n > 0) val[t] <- rnorm(n=n, mean=0, sd = sigma_theta[i])
      return(val)
    }) %>% do.call(cbind, .)
  }



  #Generate L
  sigma_L <- sqrt(1/(pi_L*J))
  if(sporadic_pleiotropy){
    L_mat <- purrr::map(seq(K), function(i){
      l <- rbinom(n=J, size=1, prob = pi_L[i])
      n <- sum(l==1)
      l[l==1] <- rnorm(n=n, mean=0, sd = sigma_L[i])
      return(l)
    }) %>% do.call(cbind, .)
  }else{
    if(sum(pi_L) > 1){
      stop("You have requested too many factors and too many causal variants to use sporadic_pleiotropy = FALSE.\n")
    }
    p <- sum(pi_L)
    nz_L_ix <- sample(c(0, seq(K)), size = J, replace = TRUE, prob = c(1-p, pi_L) )
    L_mat <- purrr::map(seq(K), function(i){
      l <- which(nz_L_ix == i)
      n <- length(l)
      val <- rep(0, J)
      if(n > 0) val[l] <- rnorm(n=n, mean=0, sd = sigma_L[i])
      return(val)
    }) %>% do.call(cbind, .)
  }


  # Compute Beta, standardized effects
  # Since phenos are scaled to variance 1, sqrt(N_m)*beta_{j,m} = z_{j,m}
  beta_std = L_mat %*% t(F_mat) + theta
  Z <- beta_std %*% diag(sqrt(N), nrow = M)

  #Compute row covariance
  if(overlap_prop > 0){
    Sigma_G <- F_mat %*% t(F_mat) + J*diag(pi_theta*sigma_theta^2)

    sigma_2_F <- (1-h_2_factor)/(h_2_factor)
    Sigma_FE <- F_mat %*% diag(sigma_2_F) %*% t(F_mat)

    if(any(h_2_trait + diag(Sigma_FE) > 1)){
      stop("Provided parameters are incompatible with generated F.\n")
    }

    sigma_E <- sqrt(1 - h_2_trait - diag(Sigma_FE))
    Sigma_E <- diag(sigma_E) %*% R_E %*% diag(sigma_E)
    #correlation of z-scores
    R <- Sigma_G + Sigma_FE + Sigma_E
    R <- overlap_prop*R + (1- overlap_prop)*diag(M)

    # Covariance of normalized effects
    # Sigma <- diag(sqrt(1/N)) %*% R %*% diag(sqrt(1/N))

    # Compute proportion of environmental variance from factors
    tau <- diag(Sigma_FE)/(1-h_2_trait)
  }else{
    R <- diag(M)
    tau <- NULL
    R_E = NULL
  }
  #Generate sampling error
  E_Z <- MASS::mvrnorm(n=J, mu = rep(0, M), Sigma = R)

  #Generate summary statistics
  if(missing(R_LD) | is.null(R_LD)){
    se_beta_hat <- matrix(1/sx) %*% matrix(1/sqrt(N), nrow = 1) # J by M
    beta_hat <- (Z + E_Z)*se_beta_hat
    # Convert L and theta from standardized scale
    L_mat <- ((1/sx)*matrix(1, nrow = J, ncol = K))*L_mat
    theta <- ((1/sx)*matrix(1, nrow = J, ncol = M))*theta

    ret <- list(beta_hat =beta_hat, se_beta_hat = se_beta_hat,
                L_mat = L_mat, F_mat = F_mat, theta = theta,
                R_E = R_E, tau = tau, R=R, Z = Z)
    return(ret)
  }

  #If LD, introduce correlation

  #Figure out how much/how many replicates of supplied LD we need
  nblock <- length(R_LD)
  l <- sapply(R_LD, function(e){length(e$values)})
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

  #snp info
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


  se_beta_hat <- matrix(1/sx) %*% matrix(1/sqrt(N), nrow = 1) # J by M
  beta_hat <- Z_hat*se_beta_hat

  # Convert L and Theta to observed scale by dividing by se of genotypes
  L_mat <- L_mat_direct <-  ((1/sx)*matrix(1, nrow = J, ncol = K))*L_mat
  theta <- theta_direct <- ((1/sx)*matrix(1, nrow = J, ncol = M))*theta

  # Transform L by LD matrix
  L_mat <- L_mat*sx # S^-inv L (the N will cancel)
  L_mat <- lapply(seq_along(block_index), function(i){
    with(R_LD[[block_index[i]]],
         vectors %*% diag(values) %*% t(vectors) %*% L_mat[start_ix[i]:end_ix[i], ])
  }) %>% do.call( rbind, .)
  L_mat <- L_mat/sx

  # Transform Theta by LD matrix
  theta <- theta*sx
  theta <- lapply(seq_along(block_index), function(i){
    with(R_LD[[block_index[i]]], vectors %*% diag(values) %*% t(vectors) %*% theta[start_ix[i]:end_ix[i], ])
  }) %>% do.call( rbind, .)
  theta <- theta/sx

  ret <- list(beta_hat =beta_hat, se_beta_hat = se_beta_hat, Z = Z,
              L_mat = L_mat, F_mat = F_mat, theta = theta,
              L_mat_direct = L_mat_direct, theta_direct = theta_direct,
              R_E = R_E, tau = tau, R = R, snp_info = snp_info_full)
  return(ret)
}







