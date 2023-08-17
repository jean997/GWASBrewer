#'@title Resample Summary Statistics for Existing Simulation Object
#'@param dat Object output by \code{sim_mv}
#'@param N Sample size, scalar, vector, matrix. See \code{?sim_mv} for more details.
#'@param R_LD LD pattern (optional). See \code{?sim_mv} for more details.
#'@param af Allele frequencies (optional, allowed only if \code{R_LD} is missing). See \code{?sim_mv} for more details.
#'@details This function can be used to generate new summary statistics for an existing simulation object.
#'The new summary statistics will have the same true causal effects as the original. However, you can use a new
#'sample size and new LD if desired. This function is primarily a wrapper for \code{gen_bhat_from_b}. This function
#'differs from \code{gen_bhat_from_b} in its arguments. It can accept the original simulated data directly rather than
#'requiring the  user to extract and supply the joint effects. It will also preserve the direct and total trait effects objects in
#'the original data.
#'
#'Note about effect scalings: If the original simulation object was generated with no allele frequencies, then
#'all effect size and effect estimate related outputs will all be on the standardized scale (in units of SD change in Y per sd change in genotype).
#'If you use this object in \code{resample_sumstats} but add allele frequencies, the resulting object will be converted
#'to non-standardized effects (effects in units of SD change in Y per alternate allele).
#'
#'If the original simulation was generated with allele frequencies, then the original object will be on the non-standardized scale.
#'\code{resample_sumstats} will sample summary statistics assuming the same non-standardized effects. This means that if you use
#'different allele frequencies with \code{resample_sumstats} than used with the original, the effects in the two objects will differ
#'on the standardized scale. The rationale for this behavior is that the non-standardized effects are the "biological" effects.
#'
#'Note that changing the allele frequency and/or LD will lead to a change in hertiability and genetic covariance (\code{Sigma_G}).
#'
#'@examples
#' # Use resample_sumstats to generate new GWAS results with the same effect sizes.
#' N <- matrix(1000, nrow = 2, ncol =2)
#' G <- matrix(0, nrow = 2, ncol = 2)
#' R_E <- matrix(c(1, 0.8, 0.8, 1), nrow = 2, ncol = 2)
#' # original data
#' dat <- sim_mv(N = N, J = 20000, h2 = c(0.4, 0.3), pi = 1000/20000,
#'                G = G, R_E = R_E)
#' # data for second GWAS
#' # Since we didn't supply af or an LD pattern in the original GWAS,
#' # we have standardized effects.
#' dat_new <- resample_sumstats(dat, N = 40000)
#'@export
resample_sumstats <- function(dat,
                              N,
                              R_LD = NULL,
                              af = NULL,
                              est_s = FALSE,
                              geno_scale = c("allele", "sd"),
                              new_h2 = NULL,
                              new_trait_corr = NULL){
  if(!"sim_mv" %in% class(dat)) stop("dat must have class sim_mv (use the sim_mv function to produce dat).\n")
  if(!is.null(R_LD) & is.null(af)) stop("af is required if R_LD is provided.\n")

  M <- ncol(dat$beta_hat)

  new_dat <- dat

  # compute variance explained by genetics in new population
  # not heritability despite using compute_h2 function
  new_dat$Sigma_G <- compute_h2(b_joint = dat$beta_joint,
                        geno_scale = dat$geno_scale,
                        pheno_sd = 1,
                        R_LD = R_LD,
                        af = af)
  v_G <- diag(new_dat$Sigma_G)
  if(!all.equal(new_dat$Sigma_G, dat$Sigma_G)){
    message("Genetic variance in the new population differs from the geneti variance in the old population.\n")
  }
  if(is.null(new_h2)){
    message("Since new_h2 is not provided, I will assume that the environmental variance is the same in the old and new population.\n")
    v_E <- diag(dat$Sigma_E)
  }else{
    new_h2 <- check_scalar_or_numeric(new_h2, "new_h2", M)
    v_E <- v_G*(1-new_h2)/new_h2
  }
  new_dat$pheno_sd <- sqrt(v_G + v_E)
  if(is.null(new_trait_corr)){
    message("Since new_trait_corr is not provided, I will assume that environmental correlation is the same in the old and new population.
            Note that this could result in different overall trait correlations.\n")
    R_E <- cov2cor(dat$Sigma_E)
    new_dat$Sigma_E <- diag(sqrt(v_E)) %*% R_E %*% diag(sqrt(v_E))
    new_dat$trait_corr <- new_dat$Sigma_G + new_dat$Sigma_E
  }else{
    new_trait_corr <- check_matrix(new_trait_corr, "new_trait_corr", M, M)
    new_trait_corr <- check_psd(new_trait_corr)
    if(!all(diag(new_trait_corr) == 1)) stop("new_trait_corr should be a correlation matrix. Found diagonal entries not equal to 1.\n")
    new_dat$trait_corr <- new_trait_corr
    Sigma_tot <- diag(new_dat$pheno_sd) %*% new_trait_corr %*% diag(new_dat$pheno_sd)
    new_dat$Sigma_E <- Sigma_tot - new_dat$Sigma_G
  }
  if(!all(out_pheno_sd == dat$pheno_sd)){
    message("Note that the phenotype in the new population has a different variance from the phenotype in the old population.\n")
    message("I will keep the phenotype on same scale as the original data, so effect sizes in the old and new object ar comparable.
            If you would like to rescale the phenotype to have variance 1, use rescale_sumstats.\n")
  }

  if(!is.null(R_LD) & is.null(af)){
    stop("Please provide af to go with R_LD.\n")
  }
  if(!is.null(af)){
    return_unit <- "allele"
  }else{
    return_unit <- "sd"
  }

  if( dat$geno_scale == "sd"){
    message("Original data have effects on the per-genotype sd scale. I will assume that per-genotype sd effects are the same in the new and old populations.\n
            You may want to consider if converting to per-allele scaling using rescale_sumstats is more appropriate for your application.\n")
    if(!is.null(af)){
      new_dat <- rescale_sumstats(new_dat, output_geno_scale = "allele", output_pheno_sd = new_dat$pheno_sd, af = af)
    }
  }

  ## Get standardized effecst for gen_bhat_from_b
  new_dat_std <- rescale_sumstats(new_dat, output_geno_scale = "sd", output_pheno_sd = 1, verbose = FALSE)
  new_ss <- gen_bhat_from_b(b_joint = new_dat_std$beta_hat, N = N, trait_corr = new_dat$trait_corr,
                            R_LD = R_LD, af = af, est_s = est_s, return_geno_unit = return_unit,
                            return_pheno_sd = new_dat$pheno_sd)


  new_dat$beta_hat <- new_ss$beta_hat
  new_dat$se_beta_hat <- new_ss$se_beta_hat
  new_dat$beta_marg <- new_ss$beta_marg # this will be different if LD changed
  if(est_s) new_dat$s_estimate <- new_ss$s_estimate
  new_dat$R <- new_ss$R
  new_dat$snp_info <- new_ss$snp_info

  new_dat <- structure(new_dat, class = c("sim_mv", "list"))
  return(new_dat)
}
