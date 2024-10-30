#'@title Resample Summary Statistics for Existing Simulation Object
#'@param dat Object output by \code{sim_mv}
#'@param N Sample size, scalar, vector, matrix. See \code{?sim_mv} for more details.
#'@param R_LD LD pattern (optional). See \code{?sim_mv} for more details.
#'@param af Allele frequencies. See \code{?sim_mv} for more details.
#'@param est_s Logical, should estimates of se(beta_hat) be produced.
#'@param geno_scale Either "allele" or "sd". Specifies the scale of the effect sizes in the output data.
#'@param new_env_var Optional. The environmental variance in the new population.
#'If missing the function will assume the environmental variance is the same as in the old population.
#'@param new_h2 Optional. The heritability in the new population. Provide at most one of \code{new_env_var} and \code{new_h2}.
#'@param new_R_E Optional, specify environmental correlation in the new population.
#'If missing, the function will assume the environmental correlation is the same as in the original data.
#'@param new_R_obs Optional, specify overall trait correlation in the new population. Specify at most one of \code{new_R_E} or \code{new_R_obs}.
#'If missing, the function will assume the environmental correlation is the same as in the original data.
#'@details This function can be used to generate new summary statistics for an existing simulation object.
#' For a discussion of this function and \code{resample_inddata}, see the "Resampling" vignette.
#'@examples
#' # Use resample_sumstats to generate new GWAS results with the same effect sizes.
#' N <- matrix(1000, nrow = 2, ncol =2)
#' G <- matrix(0, nrow = 2, ncol = 2)
#' R_E <- matrix(c(1, 0.8, 0.8, 1), nrow = 2, ncol = 2)
#' # original data
#' dat <- sim_mv(N = N, J = 20000, h2 = c(0.4, 0.3), pi = 1000/20000,
#'                G = G, R_E = R_E)
#' # data for second GWAS
#' dat_new <- resample_sumstats(dat,
#'                              N = 40000)
#'@export
resample_sumstats <- function(dat,
                              N,
                              R_LD = NULL,
                              af = NULL,
                              est_s = FALSE,
                              geno_scale = NULL,
                              new_env_var = NULL,
                              new_h2 = NULL,
                              new_R_E = NULL,
                              new_R_obs = NULL){
  if(!"sim_mv" %in% class(dat)) stop("dat must have class sim_mv (use the sim_mv function to produce dat).\n")
  if(!is.null(R_LD) & is.null(af)) stop("af is required if R_LD is provided.\n")
  if(!is.null(new_R_obs) & !is.null(new_R_E)) stop("Provide only one of new_R_obs and new_R_E.\n")
  if(!is.null(new_h2) & !is.null(new_env_var)) stop("Provide only one of new_h2 and new_env_var.\n")
  if(is.null(geno_scale)){
    geno_scale <- dat$geno_scale
  }else{
    geno_scale <- match.arg(geno_scale, choices = c("allele", "sd"))
  }
  if(geno_scale == "allele" & is.null(af)){
    stop("If geno_scale = allele, af must be supplied.\n")
  }
  M <- ncol(dat$beta_hat)
  J <- nrow(dat$beta_hat)

  new_dat <- dat

  # compute variance explained by genetics in new population
  # not heritability despite using compute_h2 function
  new_dat$Sigma_G <- compute_h2(b_joint = dat$beta_joint,
                        geno_scale = dat$geno_scale,
                        pheno_sd = 1,
                        R_LD = R_LD,
                        af = af,
                        full_mat = TRUE)
  v_G <- diag(new_dat$Sigma_G)
  if(!all(new_dat$Sigma_G == dat$Sigma_G)){
    message("Genetic variance in the new population differs from the genetic variance in the old population.")
  }
  if(is.null(new_env_var) & is.null(new_h2)){
    message("I will assume that the environmental variance is the same in the old and new population.")
    v_E <- diag(dat$Sigma_E)
  }else if(!is.null(new_env_var)){
    v_E <- check_scalar_or_numeric(new_env_var, "new_env_var", M)
  }else{
    new_h2 <- check_scalar_or_numeric(new_h2, "new_h2", M)
    v_E <- v_G*(1-new_h2)/new_h2
  }
  new_dat$pheno_sd <- sqrt(v_G + v_E)
  new_dat$h2 <- v_G/(v_G + v_E)
  if(is.null(new_R_obs) & is.null(new_R_E)){
    message("I will assume that environmental correlation is the same in the old and new population. Note that this could result in different overall trait correlations.")
    R_E <- stats::cov2cor(dat$Sigma_E)
    new_dat$Sigma_E <- diag(sqrt(v_E),nrow = M) %*% R_E %*% diag(sqrt(v_E), nrow = M)
    new_dat$trait_corr <- stats::cov2cor(new_dat$Sigma_G + new_dat$Sigma_E)
  }else if(!is.null(new_R_obs)){
    new_R_obs <- check_matrix(new_R_obs, "new_R_obs", M, M)
    new_R_obs <- check_psd(new_R_obs, "new_R_obs")
    if(!all.equal(diag(new_R_obs), rep(1, M))) stop("new_R_obs should be a correlation matrix. Found diagonal entries not equal to 1.")
    new_dat$trait_corr <- new_R_obs
    Sigma_tot <- diag(new_dat$pheno_sd, nrow  = M) %*% new_R_obs %*% diag(new_dat$pheno_sd, nrow = M)
    new_dat$Sigma_E <- Sigma_tot - new_dat$Sigma_G
    new_dat$Sigma_E <- tryCatch(check_psd(new_dat$Sigma_E, "Sigma_E"), error = function(e){
      stop("new_R_obs is incompatible with trait relationships and heritability.")
    })
  }else if(!is.null(new_R_E)){
    new_R_E <- check_matrix(new_R_E, "new_R_E", M, M)
    new_R_E <- check_psd(new_R_E, "new_R_E")
    if(!all.equal(diag(new_R_E), rep(1, M))) stop("new_R_E should be a correlation matrix. Found diagonal entries not equal to 1.")
    new_dat$Sigma_E <- diag(sqrt(v_E),nrow = M) %*% new_R_E %*% diag(sqrt(v_E), nrow = M)
    new_dat$trait_corr <- stats::cov2cor(new_dat$Sigma_G + new_dat$Sigma_E)
  }
  if(!all(new_dat$pheno_sd == dat$pheno_sd)){
    message("Note that the phenotype in the new population has a different variance from the phenotype in the old population.")
    message("I will keep the phenotype on the same scale as the original data, so effect sizes in the old and new object are comparable. If you would like to rescale the phenotype to have variance 1, use rescale_sumstats.")
  }

  if(!is.null(R_LD) & is.null(af)){
    stop("Please provide af to go with R_LD.\n")
  }
  if( dat$geno_scale == "sd"){
    message("Original data have effects on the per-genotype sd scale. I will assume that per-genotype sd effects are the same in the new and old populations.")
    if(geno_scale == "allele"){
      message("New data will be converted to the per-allele scale.")
      nr <- ceiling(J/length(af))
      new_dat <- rescale_sumstats(new_dat, output_geno_scale = "allele", output_pheno_sd = new_dat$pheno_sd, af = rep(af, nr)[1:J])
    }
  }else if(geno_scale == "sd"){
    message("New data will be converted to the per-genotype sd scale.")
    new_dat <- rescale_sumstats(new_dat, output_geno_scale = "sd", output_pheno_sd = new_dat$pheno_sd)
  }

  #new_dat_std <- rescale_sumstats(new_dat, output_geno_scale = "sd", output_pheno_sd = 1)
  ## Get standardized effects for gen_bhat_from_b
  new_ss <- gen_bhat_from_b(b_joint = new_dat$beta_joint,
                            N = N,
                            trait_corr = new_dat$trait_corr,
                            R_LD = R_LD,
                            af = af,
                            est_s = est_s,
                            input_geno_scale = new_dat$geno_scale,
                            input_pheno_sd = new_dat$pheno_sd,
                            output_geno_scale = geno_scale,
                            output_pheno_sd = new_dat$pheno_sd)



  new_dat$beta_hat <- new_ss$beta_hat
  new_dat$se_beta_hat <- new_ss$se_beta_hat
  new_dat$beta_marg <- new_ss$beta_marg # this will be different if LD changed
  if(est_s) new_dat$s_estimate <- new_ss$s_estimate
  new_dat$R <- new_ss$R
  new_dat$snp_info <- new_ss$snp_info
  new_dat$geno_scale <- geno_scale

  if(!is.null(R_LD)){
    direct_marg <- compute_R_times_mat(R_LD = R_LD, af = af, J = J, X = new_dat$direct_SNP_effects_joint*new_ss$sx)
    new_dat$direct_SNP_effects_marg <- direct_marg/new_ss$sx
    new_dat <- structure(new_dat, class = c("sim_mv", "list"))
  }
  return(new_dat)
}
