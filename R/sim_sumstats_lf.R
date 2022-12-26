
#See simulating_data.Rmd for explanation
#'@title Simulate summary statistics
#'@description Simulate summary statistics from a specified factor structure
#'@param F_mat factor matrix M by K (M = number of traits, K = number of factors)
#'@param N GWAS sample size. N can be a scalar, vector, or matrix. If N is a scalar, all GWAS have the same sample size
#'and there is no overlap between studies. If N is a vector, each element of N specifies the sample size of the corresponding
#'GWAS and there is no overlap between studies. If N is a matrix, N_ii specifies the sample size of study i
#' and N_ij specifies the number of samples present in both study i and study j. The elements of N must be positive
#' but non-integer values will not generate an error.
#'@param J Total number of SNPs to generate
#'@param h2_trait Heritability of each trait. Length M vector.
#'@param omega Proportion of trait heritability mediated by factors. Length M vector.
#'@param h2_factor Heritability of each factor. Length K vector.
#'@param pi_L Proportion of non-zero elements in L_k. Length K factor
#'@param pi_theta Proportion of non-zero elements in theta. Scalar or length M vector.
#'@param R_E Correlation between environmental trait components not mediated by factors. M by M pd matrix.
#'@param af Allele frequency (optional). A scalar, a vector, or a function that takes a number n and returns n values between 0 and 1.
#'This argument is ignored if R_LD and snp_info are provided. If af, R_LD, and snp_info
#'are all missing, SNPs are assumed to be scaled to variance 1,
#'@param R_LD List of eigen-decompositions of LD correlation matrices, may be missing.
#'@param snp_info If R_LD is provided, provide a data frame with columns "SNP" and "AF"
#'@param sporadic_pleiotropy Allow a single SNP to affect multiple factors (default TRUE).
#'@details
#'This function will generate GWAS summary statistics for M traits with K common factors.
#'The matrix F_mat provides the effects of each factor on each trait, \code{F_mat[i,j]}
#'gives the effect of factor j on trait i. The rows of \code{F_mat} will be scaled in order
#'to provide desired proportion of hertiability of each trait explained by factors but the relative
#'size and sign of elements within rows will be retained.
#'
#'Previously, if F_mat was not provided, a random factor matrix would be generated.
#'In the current version, this step must be done manually, see examples.
#'
#'With these parameters, it is possible to supply a non-feasible set of parameters.
#'Usually this occurs if the heritability of the factors is low but the heritability
#'of the traits is high leading to a contradiction. The function will return an error if this happens.
#'
#'Trait covariance: Each trait is composed of four independent components, the genetic component mediated by factors,
#'the environmental component mediated by factors, the genetic component not mediated by factors,
#'and the environmental component not mediated by factors. Therefore, the total trait covariance can be decomposed into
#'the sum of four corresponding covariance matrices. The matrix R_E specifies the correlation of the last component only.
#'In the returned object, R gives the overall trait correlation matrix multiplied by the overlap proportion matrix,
#'which is equal to the correlation in the error terms of \code{beta_hat} (see examples).
#'
#'@examples
#'myF <- generate_random_F(K = 3, M = 10, nz_factor = c(2, 3, 2),
#'                         omega = rep(0.8, 10), h2_trait = rep(0.6, 10), pad = TRUE)
#'dat <- sim_sumstats_lf(myF, N = 10000, J = 20000, h2_trait = rep(0.6, 10),
#'                       omega = rep(0.8, 10), pi_L = 0.1, pi_theta = 0.1)
#'
#'myF <- diag(2)
#'N <- matrix(c(10000, 8000, 8000, 10000), nrow = 2)
#'R_E <- matrix(c(1, 0.6, 0.6, 1), nrow = 2)
#'dat <- sim_sumstats_lf(F_mat = myF, N = N, J = 20000, h2_trait = rep(0.6, 2), omega = rep(1, 2), h2_factor = rep(1, 2),
#'                        pi_L = 0.1, pi_theta = 0.1, R_E = R_E)
#'dat$R
#'cor(dat$beta_hat[,1]-dat$beta_joint[,1], dat$beta_hat[,2]-dat$beta_joint[,2])


#'@export
sim_sumstats_lf <- function(F_mat, N, J, h2_trait, omega, h2_factor,
                            pi_L, pi_theta, R_E=NULL,
                            af = NULL, R_LD = NULL, snp_info = NULL,
                            sporadic_pleiotropy = TRUE,
                            estimate_s = FALSE){

  # Argument checks
  F_mat <- check_matrix(F_mat, "F_mat")
  M <- nrow(F_mat)
  K <- ncol(F_mat)

  h2_trait <- check_scalar_or_numeric(h2_trait, "h2_trait", M)
  h2_trait <- check_01(h2_trait)


  omega <- check_scalar_or_numeric(omega, "omega", M)
  omega <- check_01(omega)
  pi_L <- check_scalar_or_numeric(pi_L, "pi_L", K)
  pi_L <- check_01(pi_L)
  pi_theta <- check_scalar_or_numeric(pi_theta, "pi_theta", M)
  pi_theta <- check_01(pi_theta)

  if(any(omega < 1 & pi_theta == 0)){
    i <- which(omega < 1 & pi_theta == 0)
    stop(paste0("Any traits with omega < 1 (not all heritability explained by factors) must have pi_theta > 0. CHeck traits ",
                paste0(i, collapse = ",")))
  }

  #R_E
  nn <- check_N(N, M)

  if(is.null(R_E)){
    if(!Matrix::isDiagonal(nn$Nc)){
      message("R_E not provided but overlapping samples are specified. Using R_E = diag(ntrait) for no environmental covariance.")
    }
    R_E <- diag(M)
  }
  if(!missing(h2_factor)){
    h2_factor <- check_scalar_or_numeric(h2_factor, "h2_factor", K)
    h2_factor <- check_01(h2_factor)
  }else{
    if(!Matrix::isDiagonal(nn$Nc)){
      warning("h2_factors is not supplied and there are overlapping samples. Using h2_factor = 1\n")
    }
    h2_factor <- rep(1, K)
  }

  R_E <- check_matrix(R_E, "R_E", M, M)
  R_E <- check_psd(R_E, "R_E")


  #af
  if(is.null(R_LD)){
    if(is.null(af)){
      sx <- rep(1, J)
    }else if(class(af) == "function"){
      myaf <- af(J)
      sx <- sqrt(2*myaf*(1-myaf))
      af <- myaf
    }else{
      af <- check_scalar_or_numeric(af, "af", J)
      af <- check_01(af)
      sx <- sqrt(2*af*(1-af))
    }
  }


  #Re-scale F or generate it if it is missing
  srs <- omega*h2_trait # target rowSums(F^2)
  if(any(rowSums(F_mat == 0) == K & srs >0)){
      stop("One row of F is zero but corresponds to non-zero omega and h2_trait\n")
  }
  #Re scale rows of F
  scale <- sqrt(srs/rowSums(F_mat^2))
  scale[srs == 0] <- 0
  F_mat <- F_mat*scale



  #Generate theta
  sigma_theta <- sqrt( (1/(pi_theta*J)) * (1-omega)*h2_trait)
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
      stop("You have requested too many traits and too many causal variants to use sporadic_pleiotropy = FALSE.\n")
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


  # Compute standardized effects
  # Since phenos are scaled to variance 1, sqrt(N_m)*beta_{j,m} = z_{j,m}
  beta_std = L_mat %*% t(F_mat) + theta

  #Compute row covariance
  Sigma_G <- F_mat %*% t(F_mat) + J*diag(pi_theta*sigma_theta^2)

  sigma2_F <- (1-h2_factor)/(h2_factor)
  Sigma_FE <- F_mat %*% diag(sigma2_F) %*% t(F_mat)

  if(any(h2_trait + diag(Sigma_FE) > 1)){
    stop("Provided parameters are incompatible with F_mat.\n")
  }

  sigma_E <- sqrt(1 - h2_trait - diag(Sigma_FE))
  Sigma_E <- diag(sigma_E) %*% R_E %*% diag(sigma_E)
  # Trait correlation
  trait_corr <- Sigma_G + Sigma_FE + Sigma_E


  sum_stats <- gen_bhat_from_b(b_joint_std = beta_std,
                               trait_corr = trait_corr,
                               N = N,
                               R_LD = R_LD,
                               snp_info = snp_info,
                               af = af,
                               estimate_s = estimate_s,
                               L_mat_joint_std = L_mat,
                               theta_joint_std = theta)

  beta_marg <- with(sum_stats, Z*se_beta_hat)
  ret <- list(beta_hat =sum_stats$beta_hat,
              se_beta_hat = sum_stats$se_beta_hat,
              L_mat = sum_stats$L_mat,
              F_mat = F_mat,
              theta = sum_stats$theta,
              L_mat_joint = L_mat/sum_stats$sx,
              theta_joint = theta/sum_stats$sx,
              beta_joint = beta_std/sum_stats$sx,
              beta_marg = beta_marg,
              R_E = R_E,
              trait_corr = trait_corr,
              R=sum_stats$R,
              true_h2 = sum_stats$true_h2,
              af = af)
  if(estimate_s){
    ret$s_estimate <- sum_stats$s_estimate
  }

  if(!is.null(R_LD)){
    ret$snp_info <- sum_stats$snp_info
  }
  return(ret)
}







