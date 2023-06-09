
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
#'@param R_obs Observational correlation between traits. M by M pd matrix. At most one of R_E and R_obs can be specified.
#'#'@param af Optional vector of allele frequencies. If R_LD is not supplied, af can be a scalar, vector or function.
#'If af is a function it should take a single argument (n) and return a vector of n allele frequencies (See Examples).
#'If R_LD is supplied, af must be a vector with length equal to the size of the supplied LD pattern (See Examples).
#'@param sporadic_pleiotropy Allow sporadic pleiotropy between traits. Defaults to TRUE.
#'@param pi_exact If TRUE, the number of direct effect SNPs for each trait will be exactly equal to round(pi*J).
#'@param h2_exact If TRUE, the heritability of each trait will be exactly h2.
#'@param est_s If TRUE, return estimates of se(beta_hat).
#'@details
#'This function will generate GWAS summary statistics for M traits with K common factors.
#'The matrix F_mat provides the effects of each factor on each trait, \code{F_mat[i,j]}
#'gives the effect of factor j on trait i. The rows of \code{F_mat} will be scaled in order
#'to provide desired proportion of hertiability of each trait explained by factors but the relative
#'size and sign of elements within rows will be retained.
#'
#'A random factor matrix can be generated using \code{generate_random_F} (see Examples).
#'
#'It is possible to supply a non-feasible set of parameters.
#'Usually this occurs if the heritability of the factors is low but the heritability
#'of the traits is high leading to a contradiction. The function will return an error if this happens.
#'
#'Trait covariance: Each trait is composed of four independent components, the genetic component mediated by factors,
#'the environmental component mediated by factors, the genetic component not mediated by factors,
#'and the environmental component not mediated by factors. Therefore, the total trait covariance can be decomposed into
#'the sum of four corresponding covariance matrices.
#'
#'Cov(T) = Sigma_FG + Sigma_FE + Sigma_GDir + Sigma_EDir
#'
#'We assume that all cross-trait genetic sharing is explained by the factors so that Sigma_GDir is diagonal.
#'Each factor is a sum of a genetic component and an environmental components and factors are independent
#'(both genetic and environmental components) are independent across factors. This means that
#'Sigma_FG = F S_{FG} F^T and Sigma_FE = F S_{FE} F^T where S_{FG} and S_{FE} are diagonal matrices. The parameter R_E specifies
#'the correlation of the residual environmental component (i.e. R_E = cov2cor(Sigma_{EDir}).
#'Alternatively, if \code{R_obs} is specified, Sigma_EDir will be chosen to give the desired observational correlation.

#'In the returned object, \code{Sigma_G} is equal to the sum of the two genetic covariance components and \code{Sigma_E}
#'is equal to the sum of the two environmental components.
#'\code{R} gives the overall trait correlation matrix multiplied by the overlap proportion matrix,
#'which is equal to the correlation in the error terms of \code{beta_hat} (See Examples).
#'
#'@examples
#'myF <- generate_random_F(K = 3, M = 10, nz_factor = c(2, 3, 2),
#'                         omega = rep(0.8, 10), h2_trait = rep(0.6, 10), pad = TRUE)
#'dat <- sim_lf(myF, N = 10000, J = 20000, h2_trait = rep(0.6, 10),
#'                       omega = rep(0.8, 10), pi_L = 0.1, pi_theta = 0.1)
#'
#'myF <- diag(2)
#'N <- matrix(c(10000, 8000, 8000, 10000), nrow = 2)
#'R_E <- matrix(c(1, 0.6, 0.6, 1), nrow = 2)
#'dat <- sim_lf(F_mat = myF, N = N, J = 20000, h2_trait = rep(0.6, 2), omega = rep(1, 2), h2_factor = rep(1, 2),
#'                        pi_L = 0.1, pi_theta = 0.1, R_E = R_E)
#'dat$R
#'cor(dat$beta_hat[,1]-dat$beta_joint[,1], dat$beta_hat[,2]-dat$beta_joint[,2])
#'@export
sim_lf <- function(F_mat, N, J, h2_trait, omega, h2_factor,
                   pi_L, pi_theta,
                   R_E=NULL, R_obs = NULL,
                   af = NULL, R_LD = NULL,
                   sporadic_pleiotropy = TRUE,
                   est_s = FALSE,
                   snp_effect_function = "normal",
                   h2_exact = FALSE,
                   pi_exact = FALSE){

  # Argument checks
  F_mat <- check_matrix(F_mat, "F_mat")
  M <- nrow(F_mat)
  K <- ncol(F_mat)

  h2_trait <- check_scalar_or_numeric(h2_trait, "h2_trait", M)
  h2_trait <- check_01(h2_trait)


  omega <- check_scalar_or_numeric(omega, "omega", M)
  omega <- check_01(omega)
  pi_L <- check_pi(pi_L, J, K)
  pi_theta <- check_pi(pi_theta, J, M)



  f <- check_effect_function_list(snp_effect_function, M)


  if(any(omega < 1 & pi_theta == 0)){
    i <- which(omega < 1 & pi_theta == 0)
    stop(paste0("Any traits with omega < 1 (not all heritability explained by factors) must have pi_theta > 0. CHeck traits ",
                paste0(i, collapse = ",")))
  }

  #R_E
  nn <- check_N(N, M)

  if(is.null(R_E) & is.null(R_obs)){
    if(!Matrix::isDiagonal(nn$Nc)){
      message("R_E not provided but overlapping samples are specified. Using R_E = diag(ntrait) for no environmental covariance.")
    }
    R_E <- diag(M)
  }else if(!is.null(R_E)){
    if(!is.null(R_obs)){
      stop("Please provide only one of R_E and R_obs.\n")
    }
    R_E <- check_matrix(R_E, "R_E", M, M)
    R_E <- check_psd(R_E, "R_E")
    if(!all(diag(R_E) == 1)){
      stop("R_E should be a correlation matrix. All diagonal entries should be 1.")
    }
  }else if(!is.null(R_obs)){
    R_obs <- check_matrix(R_obs, "R_obs", M, M)
    R_obs <- check_psd(R_obs, "R_obs")
    if(!all(diag(R_obs) == 1)){
      stop("R_obs should be a correlation matrix. All diagonal entries should be 1.")
    }
  }
  # if(is.null(R_FE)){
  #   if(!Matrix::isDiagonal(nn$Nc)){
  #     message("R_FE not provided but overlapping samples are specified. Using R_FE = diag(nfactors) for no environmental covariance.")
  #   }
  #   R_FE <- diag(K)
  # }
  if(!missing(h2_factor)){
    h2_factor <- check_scalar_or_numeric(h2_factor, "h2_factor", K)
    h2_factor <- check_01(h2_factor)
  }else{
    if(!Matrix::isDiagonal(nn$Nc)){
      warning("h2_factors is not supplied and there are overlapping samples. Using h2_factor = 1\n")
    }
    h2_factor <- rep(1, K)
  }

  # R_FE <- check_matrix(R_FE, "R_FE", K, K)
  # R_FE <- check_psd(R_FE, "R_FE")


  #af
  if(is.null(R_LD)){
    AF <- check_af(af, J)
  }else{
    # We need af to feed to snp effect function
    if(is.null(af)) stop("Please provide af to go with R_LD.")
    l <- check_R_LD(R_LD, "l")
    AF <- check_af(af, sum(l), function_ok = FALSE)
    block_info <- assign_ld_blocks(l, J)
    AF <- AF[block_info$index]
  }


  #Re-scale rows of F
  # Scaling is such that variance of genetic component of each factor is 1
  # and total variance of each trait is 1
  srs <- omega*h2_trait # target rowSums(F^2)
  if(any(rowSums(F_mat == 0) == K & srs >0)){
      stop("One row of F is zero but corresponds to non-zero omega and h2_trait\n")
  }
  #Re scale rows of F
  scale <- sqrt(srs/rowSums(F_mat^2))
  scale[srs == 0] <- 0
  F_mat <- F_mat*scale



  #Generate theta
  theta <- sample_effects_matrix(J = J, M= M,
                                 pi = pi_theta,
                                 sigma = sqrt((1-omega)*h2_trait),
                                 f = f,
                                 sporadic_pleiotropy = sporadic_pleiotropy ,
                                 pi_exact = pi_exact,
                                 h2_exact = h2_exact,
                                 R_LD = R_LD,
                                 af = af)



  #Generate L
  L_mat <- sample_effects_matrix(J = J, M = K,
                                 pi = pi_L,
                                 sigma = 1,
                                 f = f,
                                 sporadic_pleiotropy = sporadic_pleiotropy ,
                                 pi_exact = pi_exact,
                                 h2_exact = h2_exact,
                                 R_LD = R_LD,
                                 af = af)


  # Compute standardized effects
  # Since phenos are scaled to variance 1, sqrt(N_m)*beta_{j,m} = z_{j,m}
  beta_std = L_mat %*% t(F_mat) + theta

  #Compute row (trait) covariance
  # genetic trait covariance
  #Sigma_G <- F_mat %*% t(F_mat) + diag((1-omega)*h2_trait, nrow = M)
  # Now using true genetic covariance
  Sigma_G <- compute_h2(b_joint_std = beta_std,
                        R_LD = R_LD,
                        af = af,
                        full_mat = TRUE)


  # trait covariance due to non-genetic factor component
  # sigma2_F is the variance of the environmental component of each factor
  varG_realized <- compute_h2(b_joint_std = L_mat, R_LD = R_LD, af = af, full_mat = FALSE)
  #cat(varG_realized, "\n")
  sigma2_FE <- varG_realized*(1-h2_factor)/(h2_factor)
  Sigma_FE <- F_mat %*% diag(sigma2_FE, nrow = K) %*% t(F_mat)
  #Sigma_FE <- F_mat %*% diag(sigma_FE^2, nrow = K) %*% t(F_mat)

  if(any(diag(Sigma_G) + diag(Sigma_FE) > 1)){
    cat(diag(Sigma_G), "\n")
    #cat(sigma_FE, "\n")
    cat(diag(Sigma_FE), "\n")
    stop("Provided parameters are incompatible with F_mat.\n")
  }

  # Trait covariance due to environmental contribution not mediated by factors
  if(!is.null(R_obs)){
    Sigma_E <- R_obs - Sigma_G - Sigma_FE
    Sigma_E <- try(check_psd(Sigma_E, "Sigma_E"), silent = TRUE)
    if( "try-error" %in% class(Sigma_E)){
      stop("R_obs is incompatible with F_mat and heritability.")
    }
  }else{
    sigma_E <- sqrt(1 - diag(Sigma_G) - diag(Sigma_FE))
    Sigma_E <- diag(sigma_E, nrow = M) %*% R_E %*% diag(sigma_E, nrow = M)
    #cat(diag(Sigma_E), "\n")
  }
  # Trait correlation
  trait_corr <- Sigma_G + Sigma_FE + Sigma_E


  sum_stats <- gen_bhat_from_b(b_joint_std = beta_std,
                               trait_corr = trait_corr,
                               N = N,
                               R_LD = R_LD,
                               af = af,
                               est_s = est_s,
                               L_mat_joint_std = L_mat,
                               theta_joint_std = theta)

  ret <- list(beta_hat =sum_stats$beta_hat,
              se_beta_hat = sum_stats$se_beta_hat,
              L_mat = sum_stats$L_mat,
              F_mat = F_mat,
              theta = sum_stats$theta,
              L_mat_joint = L_mat/sum_stats$sx,
              theta_joint = theta/sum_stats$sx,
              beta_joint = beta_std/sum_stats$sx,
              beta_marg =  with(sum_stats, Z*se_beta_hat),
              Sigma_G = Sigma_G,
              Sigma_E = Sigma_FE + Sigma_E,
              trait_corr = trait_corr,
              R=sum_stats$R,
              snp_info = sum_stats$snp_info)
  if(est_s){
    ret$s_estimate <- sum_stats$s_estimate
  }
  return(ret)
}
