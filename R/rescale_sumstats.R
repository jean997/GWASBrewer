#'@title Re-Scale Effects of a Simulation Object
#'@param dat A sim_mv object
#'@param output_geno_scale Desired genotype scale of output. Either "allele" or "sd".
#'@param output_pheno_sd Desired sd of phenotype, scalar or vector with length equal to the number of traits.
#'@param af If converting from sd to allele scale, provide a vector of allele frequencies.
#'@param verbose Print messages?
#'@details
#'This function can be used to change the genotype and phenotype scaling. To check the scaling
#'of the current object, look at the \code{gneo_scale} and \code{pheno_sd} elements.
#'If the current object is already on the allele scale and you desire the output to also
#'be on the allele scale, do not supply \code{af} (doing so will generate an error). If you
#'convert an "allele" scale object to an "sd" scale object, allele frequencies will be remoed.
#'@examples
#' # generate an initial data set
#' N <- matrix(10000, nrow = 2, ncol =2)
#' G <- matrix(c(0, 0.5, 0, 0), nrow = 2, ncol = 2)
#' dat <- sim_mv(N = N,
#'               G = G,
#'               J = 20000,
#'               h2 = c(0.4, 0.3),
#'               pi = 1000/20000,
#'               af = function(n){rbeta(n, 1, 5)})
#'# check scaling
#'dat$geno_scale # "allele"
#'dat$pheno_sd # 1 1
#'
#'# rescale phenotypes and convert to per-sd scale
#'dat2 <- rescale_sumstats(dat = dat,
#'                         output_geno_scale = "sd",
#'                         output_pheno_sd = c(1.5, 0.3))
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
