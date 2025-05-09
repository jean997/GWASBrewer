
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
#'@param est_s If TRUE, return estimates of se(beta_hat).
#'@param R_E Correlation between environmental trait components not mediated by factors. M by M pd matrix.
#'@param R_obs Observational correlation between traits. M by M pd matrix. At most one of R_E and R_obs can be specified.
#' @param R_LD List of LD blocks. R_LD should have class \code{list}.
#' Each element of R_LD can be either a) a matrix, b) a sparse matrix (class \code{dsCMatrix}) or c) an eigen decomposition (class \code{eigen}).
#' All elements should be correlation matrices, meaning that they have 1 on the diagonal and are positive definite.
#' @param af Optional vector of allele frequencies. If R_LD is not supplied, af can be a scalar, vector or function.
#'If af is a function it should take a single argument (n) and return a vector of n allele frequencies (See Examples).
#'If R_LD is supplied, af must be a vector with length equal to the size of the supplied LD pattern (See Examples).
#'@param sporadic_pleiotropy Allow sporadic pleiotropy between traits. Defaults to TRUE.
#'@param h2_exact If TRUE, the heritability of each trait will be exactly h2.
#'@param pi_exact If TRUE, the number of direct effect SNPs for each trait will be exactly equal to round(pi*J).
#'@param snp_effect_function_L,snp_effect_function_theta Optional function to generate variant effects in \code{L} or \code{theta}.
#'\code{snp_effect_function_L/theta} can be a single function
#' or list of functions of length equal to the number of factors/number of traits.
#' @param snp_info Optional \code{data.frame} of variant information to be passed to variant effect functions. If \code{R_LD} is
#' specified, \code{snp_info} should have number of rows equal to the size of the supplied LD pattern. Otherwise \code{snp_info}
#' should have \code{J} rows.
#'@details
#'This function will generate GWAS summary statistics for M traits with K common factors.
#'The matrix F_mat provides the effects of each factor on each trait, \code{F_mat[i,j]}
#'gives the effect of factor j on trait i. The rows of \code{F_mat} will be scaled in order
#'to provide desired proportion of heritability of each trait explained by factors but the relative
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
#'\deqn{Cov(T) = Sigma_{FG} + Sigma_{FE} + Sigma_{GDir} + Sigma_{EDir}}
#'
#'We assume that all cross-trait genetic sharing is explained by the factors so that \eqn{Sigma_{GDir}} is diagonal.
#'Each factor is a sum of a genetic component and an environmental components and factors are independent
#'(both genetic and environmental components) are independent across factors. This means that
#'\eqn{Sigma_{FG} = F S_{FG} F^T} and \eqn{Sigma_{FE} = F S_{FE} F^T} where \eqn{S_{FG}} and \eqn{S_{FE}} are diagonal matrices. The parameter \code{R_E} specifies
#'the correlation of the residual environmental component (i.e. \eqn{R_E = cov2cor(Sigma_{EDir}}).
#'Alternatively, if \code{R_obs} is specified, \eqn{Sigma_{EDir}} will be chosen to give the desired observational correlation.

#'In the returned object, \code{Sigma_G} is equal to the sum of the two genetic covariance components and \code{Sigma_E}
#'is equal to the sum of the two environmental components.
#'\code{R} gives the overall trait correlation matrix multiplied by the overlap proportion matrix,
#'which is equal to the correlation in the error terms of \code{beta_hat} (See Examples).
#'
#'@examples
#'myF <- generate_random_F(K = 3, M = 10, nz_factor = c(2, 3, 2),
#'                         omega = rep(0.8, 10),
#'                         h2_trait = rep(0.6, 10), pad = TRUE)
#'dat <- sim_lf(myF, N = 10000, J = 20000, h2_trait = rep(0.6, 10),
#'                       omega = rep(0.8, 10), pi_L = 0.1, pi_theta = 0.1)
#'
#'myF <- diag(2)
#'N <- matrix(c(10000, 8000, 8000, 10000), nrow = 2)
#'R_E <- matrix(c(1, 0.6, 0.6, 1), nrow = 2)
#'dat <- sim_lf(F_mat = myF, N = N, J = 20000, h2_trait = rep(0.6, 2),
#'              omega = rep(1, 2), h2_factor = rep(1, 2),
#'              pi_L = 0.1, pi_theta = 0.1, R_E = R_E)
#'dat$R
#'cor(dat$beta_hat[,1]-dat$beta_joint[,1], dat$beta_hat[,2]-dat$beta_joint[,2])
#'@export
sim_lf <- function(F_mat,
                   N,
                   J,
                   h2_trait,
                   omega,
                   h2_factor,
                   pi_L,
                   pi_theta,
                   est_s = FALSE,
                   R_E=NULL,
                   R_obs = NULL,
                   R_LD = NULL,
                   af = NULL,
                   sporadic_pleiotropy = TRUE,
                   h2_exact = FALSE,
                   pi_exact = FALSE,
                   snp_effect_function_L = "normal",
                   snp_effect_function_theta = "normal",
                   snp_info = NULL){

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

  if(any(omega < 1 & pi_theta == 0)){
    i <- which(omega < 1 & pi_theta == 0)
    stop(paste0("Any traits with omega < 1 (not all heritability explained by factors) must have pi_theta > 0. Check traits ",
                paste0(i, collapse = ",")))
  }

  #R_E
  nn <- check_N(N, M)

  if(is.null(R_E) & is.null(R_obs)){
    if(nn$overlap){
      message("Neither R_E nor R_obs was provided but there are overlapping samples.\n
              I will assume that environmental components not mediated by factors are independent.")
    }
    R_E <- diag(M)
  }else if(!is.null(R_E)){
    if(!is.null(R_obs)){
      stop("Please provide only one of R_E and R_obs.\n")
    }
    R_E <- check_matrix(R_E, "R_E", M, M)
    R_E <- check_psd(R_E, "R_E")
    if(!all.equal(diag(R_E), rep(1, M), check.names = FALSE)){
      stop("R_E should be a correlation matrix. All diagonal entries should be 1.")
    }
  }else if(!is.null(R_obs)){
    R_obs <- check_matrix(R_obs, "R_obs", M, M)
    R_obs <- check_psd(R_obs, "R_obs")
    if(!all.equal(diag(R_obs), rep(1, M), check.names = FALSE)){
      stop("R_obs should be a correlation matrix. All diagonal entries should be 1.")
    }
  }

  if(!missing(h2_factor)){
    h2_factor <- check_scalar_or_numeric(h2_factor, "h2_factor", K)
    h2_factor <- check_01(h2_factor, "h2_factor")
  }else{
    if(nn$overlap){
      warning("h2_factor is not supplied and there are overlapping samples.
              Using h2_factor = 1 (no environmental variance mediated by factors).\n")
    }
    h2_factor <- rep(1, K)
  }

  #af and snp_info, these go into effect size function and we will need them
  # later for LD
  if(is.null(R_LD)){
    af <-  check_af(af, J)
    if(is.null(snp_info) & is.null(af)){
      snp_info <- data.frame(SNP = 1:J, AF = NA)
    }else if(is.null(snp_info)){
      snp_info <- data.frame(SNP = 1:J, AF = af)
    }else{
      snp_info <- check_snpinfo(snp_info, J)
      snp_info$SNP <- 1:J
      if(is.null(af)) snp_info$AF <- NA
        else snp_info$AF <- af
    }
  }else{
    if(is.null(af)) stop("Please provide af to go with R_LD.")
    l <- check_R_LD(R_LD, "l") # vector of LD block sizes
    af <- check_af(af, sum(l), function_ok = FALSE)
    block_info <- assign_ld_blocks(l, J)

    if(is.null(snp_info)){
      snp_info <- data.frame(SNP = 1:sum(l), AF = af)
    }else{
      snp_info <- check_snpinfo(snp_info, sum(l))
      snp_info$SNP <- 1:sum(l)
      snp_info$AF <- af
    }

    snp_info$block <- rep(seq(length(l)), l)
    snp_info$ix_in_block <- sapply(l, function(nl){seq(nl)}) %>% unlist()
    snp_info <- snp_info[block_info$index,] # now snp_info has J rows
    snp_info$rep <- rep(block_info$rep, block_info$l)
    snp_info$SNP <- with(snp_info, paste0(SNP, ".", rep))
  }

  f_L <- check_effect_function_list(snp_effect_function_L, K, snp_info)
  f_theta <- check_effect_function_list(snp_effect_function_theta, M, snp_info)


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
                                 f = f_theta,
                                 sporadic_pleiotropy = sporadic_pleiotropy ,
                                 pi_exact = pi_exact,
                                 h2_exact = h2_exact,
                                 snp_info = snp_info,
                                 R_LD = R_LD,
                                 af = af)



  #Generate L
  #message(paste0("h2_trait: ", paste0(h2_trait, collapse = ",")))
  L_mat <- sample_effects_matrix(J = J, M = K,
                                 pi = pi_L,
                                 sigma = as.numeric(h2_factor > 0 & colSums(F_mat^2) > 0),
                                 f = f_L,
                                 sporadic_pleiotropy = sporadic_pleiotropy ,
                                 pi_exact = pi_exact,
                                 h2_exact = h2_exact,
                                 snp_info = snp_info,
                                 R_LD = R_LD,
                                 af = af)


  # Compute standardized effects
  # Since phenos are scaled to variance 1, sqrt(N_m)*beta_{j,m} = z_{j,m}
  beta_std = L_mat %*% t(F_mat) + theta

  #Compute row (trait) covariance
  # genetic trait covariance
  #Sigma_G <- F_mat %*% t(F_mat) + diag((1-omega)*h2_trait, nrow = M)
  # Now using true genetic covariance
  Sigma_G <- compute_h2(b_joint = beta_std,
                        geno_scale = "sd",
                        pheno_sd = 1,
                        R_LD = R_LD,
                        af = af,
                        full_mat = TRUE)


  # trait covariance due to non-genetic factor component
  # realized genetic variance of each factor, close to 1
  varG_realized <- compute_h2(b_joint = L_mat,
                              geno_scale = "sd",
                              pheno_sd = 1,
                              R_LD = R_LD,
                              af = af,
                              full_mat = FALSE)
  #cat(varG_realized, "\n")
  sigma2_FE <- varG_realized*(1-h2_factor)/(h2_factor) # Variance from environmental components of factors
  Sigma_FE <- F_mat %*% diag(sigma2_FE, nrow = K) %*% t(F_mat)

  if(any(diag(Sigma_G) + diag(Sigma_FE) > 1)){
    cat(diag(Sigma_G), "\n")
    #cat(sigma_FE, "\n")
    cat(diag(Sigma_FE), "\n")
    stop("Provided parameters are incompatible with F_mat.\n")
  }

  # Trait covariance due to environmental contribution not mediated by factors
  if(!is.null(R_obs)){
    Sigma_E <- R_obs - Sigma_G - Sigma_FE
    Sigma_E <- tryCatch(check_psd(Sigma_E, "Sigma_E"), error = function(e){
      stop("R_obs is incompatible with trait relationships and heritability.")
    })
  }else{
    sigma_E <- sqrt(1 - diag(Sigma_G) - diag(Sigma_FE))
    Sigma_E <- diag(sigma_E, nrow = M) %*% R_E %*% diag(sigma_E, nrow = M)
  }
  # Trait correlation
  trait_corr <- Sigma_G + Sigma_FE + Sigma_E

  sum_stats <- gen_bhat_from_b(b_joint = beta_std,
                               trait_corr = trait_corr,
                               N = N,
                               R_LD = R_LD,
                               af = af,
                               est_s = est_s,
                               input_geno_scale = "sd",
                               output_geno_scale = case_when(is.null(af) ~ "sd",
                                                            TRUE ~ "allele"),
                               input_pheno_sd = 1,
                               output_pheno_sd = 1)


  sum_stats$F_mat <- F_mat
  sum_stats$Sigma_G <- Sigma_G
  sum_stats$Sigma_E <- Sigma_FE + Sigma_E
  sum_stats$Sigma_FE <- Sigma_FE
  sum_stats$h2 <- diag(Sigma_G)
  sum_stats$trait_corr <- trait_corr
  sum_stats$snp_info <- snp_info

  if(is.null(af)){
    sum_stats$L_mat_joint <-  L_mat
    sum_stats$theta_joint <- theta
  }else{
    sum_stats$L_mat_joint <- L_mat/sum_stats$sx
    sum_stats$theta_joint <- theta/sum_stats$sx
  }
  if(is.null(R_LD)){
    sum_stats$L_mat_marg <- sum_stats$L_mat_joint
    sum_stats$theta_marg <- sum_stats$theta_joint
  }else{
    sum_stats$L_mat_marg <- compute_R_times_mat(R_LD, af, J, L_mat)
    sum_stats$L_mat_marg <- sum_stats$L_mat_marg/sum_stats$sx
    sum_stats$theta_marg <- compute_R_times_mat(R_LD, af, J, theta)
    sum_stats$theta_marg <- sum_stats$theta_marg/sum_stats$sx
  }
  sum_stats <- structure(sum_stats, class = c("sim_lf", "list"))

  return(sum_stats)
}
