#'@export
rescale_sumstats <- function(dat,
                             output_geno_scale = c("allele", "sd"),
                             output_pheno_sd = 1,
                             af = NULL){

  M <- ncol(dat$beta_hat)
  output_pheno_sd <- check_scalar_or_numeric(output_pheno_sd, "out_pheon_sd", M)
  if(! all(output_pheno_sd > 0)){
    stop("Phenotype SD must be greater than or equal to 0.\n")
  }
  J <- nrow(dat$beta_hat)
  if(dat$geno_scale != output_geno_scale){
    if(dat$geno_scale == "allele"){
      message("Converting effects from per-allele to per-genotype SD.\n")
      if(any(is.na(dat$snp_info$AF))){
        stop("Something is wrong. dat$geno_scale is allele but AF is missing from dat$snp_info.\n")
      }
      sx <- with(dat$snp_info, sqrt(2*AF*(1-AF)))
      dat$beta_hat <- col_times(dat$beta_hat, sx)
      dat$se_beta_hat <- col_times(dat$se_beta_hat, sx)
      dat$beta_joint <- col_times(dat$beta_joint, sx)
      dat$beta_marg <- col_times(dat$beta_marg, sx)

      dat$L_mat <- col_times(dat$L_mat, sx)
      dat$theta <- col_times(dat$theta, sx)

      dat$s_estimate <- col_times(dat$s_estimate, sx)

      dat$direct_SNP_effects_joint <- col_times(dat$direct_SNP_effects_joint, sx)
      dat$direct_SNP_effects_joint <- col_times(dat$direct_SNP_effects_joint, sx)

      dat$geno_scale <- "sd"
      dat$snp_info$AF <- NA
    }else if(dat$geno_scale == "sd"){
      if(is.null(af)){
        stop("Please provide af to convert to per-allele effect scale.\n")
      }
      message("Converting effects from per-genotype SD to per-allele.\n")
      af <- check_af(af, J)
      sx <- sqrt(2*af*(1-af))
      dat$snp_info$AF <- af
      dat$beta_hat <- col_times(dat$beta_hat, 1/sx)
      dat$se_beta_hat <- col_times(dat$se_beta_hat, 1/sx)
      dat$beta_joint <- col_times(dat$beta_joint, 1/sx)
      dat$beta_marg <- col_times(dat$beta_marg, 1/sx)

      dat$L_mat <- col_times(dat$L_mat, 1/sx)
      dat$theta <- col_times(dat$theta, 1/sx)

      dat$s_estimate <- col_times(dat$s_estimate, 1/sx)

      dat$direct_SNP_effects_joint <- col_times(dat$direct_SNP_effects_joint, 1/sx)
      dat$direct_SNP_effects_joint <- col_times(dat$direct_SNP_effects_joint, 1/sx)

      dat$geno_scale <- "allele"
    }
  }
  if(! identical(dat$pheno_sd, output_pheno_sd)){
    message("Converting effects to change scale of phenotypes.\n")
    pheno_scale <- output_pheno_sd/dat$pheno_sd
    dat$beta_hat <- row_times(dat$beta_hat, pheno_scale)
    dat$se_beta_hat <- row_times(dat$se_beta_hat, pheno_scale)
    dat$beta_joint <- row_times(dat$beta_joint, pheno_scale)
    dat$beta_marg <- row_times(dat$beta_marg, pheno_scale)

    dat$L_mat <- row_times(dat$L_mat, pheno_scale)
    dat$theta <- row_times(dat$theta, pheno_scale)

    dat$s_estimate <- row_times(dat$s_estimate, pheno_scale)

    dat$direct_SNP_effects_joint <- row_times(dat$direct_SNP_effects_joint, pheno_scale)
    dat$direct_SNP_effects_joint <- row_times(dat$direct_SNP_effects_joint, pheno_scale)


    dat$pheno_sd <- output_pheno_sd
  }
  return(dat)

}

row_times <- function(A, b){
  if(is.null(A)) return(NULL)
  t(t(A)*b)
}

col_times <- function(A, b){
  if(is.null(A)) return(NULL)
  A*b
}
