#'@title Compute heritability from standardized or non-standardized effects
#'@param b_joint_std,b_joint matrix of standardized or non-standardized effects. Provide only one of these options.
#'@param R_LD LD pattern (optional). See \code{?sim_mv} for more details.
#'@param af Allele frequencies (optional, allowed only if \code{R_LD} is missing). See \code{?sim_mv} for more details.
#'@param full_mat If TRUE, return the full genetic variance-covariance matrix
#'@export
compute_h2 <- function(b_joint,
                       geno_scale = c("allele", "sd"),
                       pheno_sd = 1,
                       R_LD = NULL,
                       af = NULL,
                       full_mat = FALSE){

  M <- ncol(b_joint)
  J <- nrow(b_joint)

  pheno_sd <- check_scalar_or_numeric(pheno_sd, "pheno_sd", M)
  if(!all(pheno_sd == 1)){
    b_joint <- row_times(b_joint, 1/pheno_sd)
  }

  if(is.null(R_LD)){
    if(geno_scale == "allele"){
      if(is.null(af)){
        stop("Provide af if geno_scale = allele.\n")
      }
      af <- check_af(af, J)
      sx <- sqrt(2*af*(1-af))
      b_joint_std <- b_joint*sx
    }else{
      b_joint_std <- b_joint
    }
    if(!full_mat) return(colSums(b_joint_std^2))
    return(t(b_joint_std) %*% b_joint_std)
  }

  ld_mat <- check_R_LD(R_LD, "matrix")
  if(is.null(af)){
    stop("Provide af if R_LD is not NULL.\n")
  }
  l <- check_R_LD(R_LD, "l")
  af_ld <- check_af(af, sum(l), function_ok = FALSE)
  nblock <- length(l)
  block_info <- assign_ld_blocks(l, J)
  if(!is.null(block_info$last_block_info)){
    b <- block_info$last_block_info[1]
    x <- block_info$last_block_info[2]
    last_block <- ld_mat[[b]][seq(x), seq(x)]
    ld_mat[[nblock + 1]] <- last_block

  }
  af <- af_ld[block_info$index]
  sx <- sqrt(2*af*(1-af))
  if(geno_scale == "allele"){
    b_joint_std <- col_times(b_joint, sx)
  }else{
    b_joint_std <- b_joint
  }
  nb <- length(block_info$l)
  start_ix <- cumsum(c(1, block_info$l[-nb]))
  end_ix <- start_ix + block_info$l-1

  if(!full_mat){
    true_h2 <- lapply(seq(nb), function(i){
      # emulator::quad.diag(R, b_joint_std)
      colSums(crossprod(ld_mat[[block_info$block_index[i]]], b_joint_std[start_ix[i]:end_ix[i], ])*b_joint_std[start_ix[i]:end_ix[i], ])
    }) %>% Reduce("+", .)

    return(true_h2)
  }
  true_h2 <- lapply(seq(nb), function(i){
    # quad.form(R, b_joint_std)
    crossprod(crossprod(ld_mat[[block_info$block_index[i]]], b_joint_std[start_ix[i]:end_ix[i], ]), b_joint_std[start_ix[i]:end_ix[i], ])
  }) %>% Reduce("+", .)

  return(true_h2)

}
