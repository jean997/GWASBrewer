#'@export
rescale_sumstats <- function(dat,
                             output_geno_scale = c("allele", "sd"),
                             output_pheno_sd = 1,
                             af = NULL,
                             verbose = TRUE){

  M <- ncol(dat$beta_hat)
  output_pheno_sd <- check_scalar_or_numeric(output_pheno_sd, "out_pheno_sd", M)
  output_geno_scale <- match.arg(output_geno_scale)
  if(! all(output_pheno_sd > 0)){
    stop("Phenotype SD must be greater than or equal to 0.\n")
  }
  J <- nrow(dat$beta_hat)

  if(dat$geno_scale == "allele"){
    if(!is.null(af)){
      stop("If dat is on the allele genotype scale, the af argument must be missing. rescale_sumstats cannot be used to change allele frequencies. Use resample_sumstats instead.\n")
    }
    af <- dat$snp_info$AF
  }else{
    if(output_geno_scale == "allele" & is.null(af)){
      stop("To convert an SD scale object to allele scale, provide af.\n")
    }else if(output_geno_scale == "sd" & !is.null(af)){
      stop("dat is on the SD scale and output_geno_scale = sd, so af will not be used. Please re-run omitting this argument.\n")
    }
  }
  af <- check_scalar_or_numeric(af, "af", J)

  A <- get_convert_matrix(input_geno_scale = dat$geno_scale,
                          output_geno_scale = output_geno_scale,
                          input_pheno_sd = dat$pheno_sd,
                          output_pheno_sd = output_pheno_sd,
                          J = J,
                          af = af)

  dat$beta_hat <- mat_times(dat$beta_hat, A)
  dat$se_beta_hat <- mat_times(dat$se_beta_hat, A)
  dat$beta_joint <- mat_times(dat$beta_joint, A)
  dat$beta_marg <- mat_times(dat$beta_marg, A)

  dat$L_mat <- mat_times(dat$L_mat, A)
  dat$theta <- mat_times(dat$theta, A)

  dat$s_estimate <- mat_times(dat$s_estimate, A)

  dat$direct_SNP_effects_joint <- mat_times(dat$direct_SNP_effects_joint, A)
  dat$direct_SNP_effects_marg <- mat_times(dat$direct_SNP_effects_marg, A)

  dat$geno_scale <- output_geno_scale
  if(is.null(af)){
    dat$snp_info$AF <- NA
  }else{
    dat$snp_info$AF <- af
  }


  if(! identical(dat$pheno_sd, output_pheno_sd)){
    # Convert genetic and environmental variance - covariance matrices
    pheno_scale <- output_pheno_sd/dat$pheno_sd
    B <- kronecker(pheno_scale, pheno_scale) |> matrix(nrow = M)
    dat$Sigma_G <- B*dat$Sigma_G
    dat$Sigma_E <- B*dat$Sigma_E
    C <- kronecker(pheno_scale, 1/pheno_scale) |> matrix(nrow = M)
    dat$direct_trait_effects <- C*dat$direct_trait_effects
    dat$total_trait_effects <- C*dat$total_trait_effects
    dat$pheno_sd <- output_pheno_sd
  }
  return(dat)

}

mat_times <- function(A, B){
  if(is.null(A)) return(NULL)
  B*A
}
row_times <- function(A, b){
  if(is.null(A)) return(NULL)
  t(t(A)*b)
}

col_times <- function(A, b){
  if(is.null(A)) return(NULL)
  A*b
}

get_convert_matrix <- function(input_geno_scale,
                               output_geno_scale,
                               input_pheno_sd,
                               output_pheno_sd,
                               J,
                               af = NULL){
  ## no checks, internal function
  M <- length(input_pheno_sd)
  if(is.null(af)){
    sx <- rep(1, J)
  }else{
    sx <- sqrt(2*af*(1-af))
  }

  if(identical(input_geno_scale, output_geno_scale)  & identical(input_pheno_sd, output_pheno_sd)){
    return(1)
  }else if(identical(input_geno_scale, output_geno_scale)){
    pheno_scale <- output_pheno_sd/input_pheno_sd
    scale_matrix <- matrix(pheno_scale, nrow = J, ncol = M, byrow = T)
  }else if(input_geno_scale == "sd"){ ## input is sd, output is allele
    pheno_scale <- output_pheno_sd/input_pheno_sd
    scale_matrix <- kronecker(1/sx, pheno_scale) |> matrix(nrow = J, byrow = T)
  }else{ ## input is allele, output is sd
    pheno_scale <- output_pheno_sd/input_pheno_sd
    scale_matrix <- kronecker(sx, pheno_scale) |> matrix(nrow = J, byrow = T)
  }
  return(scale_matrix)
}
