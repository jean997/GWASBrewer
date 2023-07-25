#'@title Resample Summary Statistics for Existing Simulation Object
#'@param dat Object output by \code{sim_mv}
#'@param N Sample size, scalar, vector, matrix. See \code{?sim_mv} for more details.
#'@param R_LD LD pattern (optional). See \code{?sim_mv} for more details.
#'@param af Allele frequencies (optional, allowed only if \code{R_LD} is missing). See \code{?sim_mv} for more details.
#'@details This function can be used to generate new summary statistics for an existing simulation object.
#'The new summary statistics will have the same true causal effects as the original. However, you can use a new
#'sample size and new LD if desired. This function is primarily a wrapper for \code{gen_bhat_from_b}. This function
#'differs from \code{gen_bhat_from_b} in its arguments. It can accept the original simulated data directly rather than
#'requiring the  user to extract and supply the joint effects. It will also preserver the direct and total trait effects objects in
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
                              est_s = FALSE){
  if(!"sim_mv" %in% class(dat)) stop("dat must have class sim_mv (use the sim_mv function to produce dat).\n")
  if(!is.null(R_LD) & is.null(af)) stop("af is required if R_LD is provided.\n")
  new_obj <- dat
  if( "sim_mv_std" %in% class(dat)){
    # dat contains standardized effects
    new_ss <- gen_bhat_from_b(b_joint_std = dat$beta_joint,
                              trait_cor = dat$trait_cor,
                              N = N,
                              R_LD = R_LD,
                              af = af,
                              est_s = est_s)
    if(!is.null(af)){
      warning("Original data were on standardized scale but resampled data will be on the non-standardized scale (see help page for more information).")
      # convert all objects not in new_ss to non-standardized scale
      new_obj$direct_SNP_effects_marg <- new_obj$direct_SNP_effects_marg/new_ss$sx
      new_obj$direct_SNP_effects_joint <- new_obj$direct_SNP_effects_joint/new_ss$sx
      new_obj$beta_joint <- new_obj$beta_joint/new_ss$sx

    }
  }else{
    new_ss <- gen_bhat_from_b(b_joint = dat$beta_joint,
                              trait_cor = dat$trait_cor,
                              N = N,
                              R_LD = R_LD,
                              af = af,
                              est_s = est_s)
    if(is.null(af)){
      warning("Origingal data were on the non-standardized scale but resampled data will be on the standardized scale because af was not provided (see help page for more information).")
      # convert all objects not in new_ss to standardized scale
      sx_orig <- with(dat$snp_info, sqrt(2*AF*(1-AF)))
      new_obj$direct_SNP_effects_marg <- new_obj$direct_SNP_effects_marg*sx_orig
      new_obj$direct_SNP_effects_joint <- new_obj$direct_SNP_effects_joint*sx_orig
      new_obj$beta_joint <- new_obj$beta_joint*sx_orig
    }
  }
  new_obj$beta_hat <- new_ss$beta_hat
  new_obj$se_beta_hat <- new_ss$se_beta_hat
  new_obj$beta_marg <- new_ss$beta_marg
  if(est_s) new_obj$s_estimate <- new_ss$s_estimate
  new_obj$R <- new_ss$R
  new_obj$snp_info <- new_ss$snp_info
  if(is.null(af)){
    new_obj <- structure(new_obj, class = c("sim_mv", "sim_mv_std", "list"))
    new_obj$Sigma_G <- compute_h2(b_joint_std = new_obj$beta_joint,
                                  full_mat = TRUE)
  }else{
    new_obj <- structure(new_obj, class = c("sim_mv", "list"))
    new_obj$Sigma_G <- compute_h2(b_joint = new_obj$beta_joint,
                                  R_LD = R_LD, af = af,
                                  full_mat = TRUE)
  }
  # recaclucate heritability/genetic covariance matrix with (possibly) new LD and AF
  new_obj$Sigma_E <- new_obj$trait_cor - new_obj$Sigma_G
  return(new_obj)
}
