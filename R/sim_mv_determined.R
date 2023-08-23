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
      my_af <- rep(af, ceiling(J/length(af)))[1:J]
      my_af <- check_scalar_or_numeric(my_af, "af", J)
    }
    if(is.null(af)){
      stop("af must be provided if using per-allele scale effects")
    }
  }else{
    my_af <- NULL
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
                            af = my_af)
    direct_SNP_effects_joint <- A*direct_SNP_effects_joint
    pheno_scale <- 1/output_pheno_sd
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
