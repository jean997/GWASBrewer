#'@title Simulate multivariate GWAS Data with Specified Direct Effects
#'@param N GWAS sample size. N can be a scalar, vector, or matrix. If N is a scalar, all GWAS have the same sample size
#'and there is no overlap between studies. If N is a vector, each element of N specifies the sample size of the corresponding
#'GWAS and there is no overlap between studies. If N is a matrix, N_ii specifies the sample size of study i
#' and N_ij specifies the number of samples present in both study i and study j. The elements of N must be positive
#' but non-integer values will not generate an error.
#'@param direct_SNP_effects_joint Matrix of direct variant effects. Should be variants by traits.
#'@param geno_scale Genotype scale of provided effects. Either "allele" or "sd".
#'@param pheno_sd Phenotype standard deviation, a scalar or vector of length number of traits.
#'@param G Matrix of direct effects. Rows correspond to the 'from' trait
#'and columns correspond to the 'to' trait, so \code{G[1,2]} is the direct effect of trait 1 on trait 2. G should have 0
#'on the diagonal. Be sure that \code{G} is on the same scale as the effect sizes.
#'@param est_s If TRUE, return estimates of se(`beta_hat`). If FALSE, the exact standard error of `beta_hat` is returned. Defaults to FALSE.
#'@param R_obs Total observational correlation between traits. R_obs won't impact summary statistics unless there is sample overlap.
#'See Details for default behavior.
#'@param R_E Total correlation of the environmental components only. R_E and R_obs are alternative methods of specifying trait correlation.
#'Use only one of these two options.
#'R_E may be phased out in the future.
#'@param R_LD Optional list of LD blocks. R_LD should have class \code{list}.
#'Each element of R_LD can be either a) a matrix, b) a sparse matrix (class \code{dsCMatrix}) or c) an eigen decomposition (class \code{eigen}).
#'All elements should be correlation matrices, meaning that they have 1 on the diagonal and are positive definite. See Details and vignettes.
#'@param af Optional vector of allele frequencies. If R_LD is not supplied, af can be a scalar, vector or function.
#'If af is a function it should take a single argument (n) and return a vector of n allele frequencies (See Examples).
#'If R_LD is supplied, af must be a vector with length equal to the size of the supplied LD pattern (See Examples).
#'@return A \code{sim_mv} function. See \code{?sim_mv} for details.
#'
#'@details A wrapper for \code{sim_mv}. See \code{?sim_mv} and the "Providing an Exact Set of Direct Effects" section of the Effect Size vignette.
#'
#'@examples
#' G <- matrix(c(0, 0.5, 0, 0), nrow = 2, byrow =TRUE)
#' my_effects <- matrix(0, nrow = 10, ncol = 2)
#' my_effects[c(1, 5),1] <- c(-0.008, 0.01)
#' my_effects[c(3, 6, 9), 2] <- c(-0.02, 0.06, 0.009)
#' my_effects
#' # for fun, lets include some sample overlap
#' N <- matrix(c(40000, 10000, 10000, 20000), nrow = 2)
#' sim_dat <- sim_mv_determined(N = N,
#'                               direct_SNP_effects_joint = my_effects,
#'                               geno_scale = "sd",
#'                               pheno_sd = 1,
#'                               G=G,
#'                               est_s = TRUE)
#'
#' sim_dat$direct_SNP_effects_joint
#' sim_dat$beta_joint
#' sim_dat$Sigma_G

#'@export
sim_mv_determined <- function(N,
                              direct_SNP_effects_joint,
                              geno_scale,
                              pheno_sd,
                              G=0,
                              est_s = FALSE,
                              R_obs = NULL,
                              R_E = NULL,
                              R_LD = NULL,
                              af = NULL){


  direct_SNP_effects_joint <- check_matrix(direct_SNP_effects_joint, "direct_SNP_effects_joint")
  M <- ncol(direct_SNP_effects_joint)
  J <- nrow(direct_SNP_effects_joint)
  nn <- check_N(N, M)
  geno_scale <- match.arg(geno_scale, choices = c("allele", "sd"))
  if(geno_scale == "allele"){
    if(is.null(R_LD)){
      af <- full_af <- check_scalar_or_numeric(af, "af", J)
    }else{
      full_af <- rep(af, ceiling(J/length(af)))[1:J]
      full_af <- check_scalar_or_numeric(full_af, "af", J)
    }
    if(is.null(af)){
      stop("af must be provided if using per-allele scale effects")
    }
  }else{
    full_af <- NULL
  }
  pheno_sd <- check_scalar_or_numeric(pheno_sd, "pheno_sd", M)
  if(!all(pheno_sd == 1) | geno_scale == "allele"){
    G <- check_matrix(G, "G", M, M)
    message("Converting both G and direct effects to standaridzed scale.")
    A <- get_convert_matrix(input_geno_scale = geno_scale,
                            output_geno_scale = "sd",
                            input_pheno_sd = pheno_sd,
                            output_pheno_sd = 1,
                            J = J,
                            af = full_af)
    direct_SNP_effects_joint <- A*direct_SNP_effects_joint
    pheno_scale <- 1/pheno_sd
    C <- kronecker(pheno_scale, 1/pheno_scale) |> matrix(nrow = M)
    G <- C*G
  }
  G <- check_G(G, rep(1, M))
  # step 1 caclucate heritability
  # 1a caclucate total effects
  beta_joint <- direct_SNP_effects_joint%*%(diag(1, nrow = M) + G$G_tot)
  Sigma_G <- compute_h2(b_joint = beta_joint,
                   geno_scale = "sd",
                   pheno_sd = 1,
                   R_LD = R_LD,
                   af = af,
                   full_mat = TRUE)
  h2 <- diag(Sigma_G)

  ## generate snp effect functions
  myfuncs <- lapply(1:M, function(i){
    fixed_to_scale_fam(direct_SNP_effects_joint[,i])
  })

  ## call sim_mv
  dat <- sim_mv(N = N,
                J = J,
                h2 = h2,
                pi = 1,
                G=G$G_dir,
                est_s = est_s,
                R_obs = R_obs,
                R_E = R_E,
                R_LD = R_LD,
                af = af,
                snp_effect_function = myfuncs,
                snp_info = NULL,
                sporadic_pleiotropy = TRUE,
                pi_exact = FALSE,
                h2_exact = FALSE)
  return(dat)
}
