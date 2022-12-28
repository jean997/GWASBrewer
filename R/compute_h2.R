compute_h2 <- function(b_joint_std,
                       b_joint,
                       R_LD = NULL,
                       af = NULL,
                       full_mat = FALSE){
  if(!missing(b_joint)){
    if(!missing(b_joint_std)){
      stop("Please provide only one of b_joint or b_joint_std")
    }
    if(is.null(af)){
      stop("af is required if using non-standardized effects")
    }
    b_joint <- check_matrix(b_joint, "b_joint")
    M <- ncol(b_joint)
    J <- nrow(b_joint)
    b_type <- "non_std"
  }else if(!missing(b_joint_std)){
    b_joint_std <- check_matrix(b_joint_std, "b_joint_std")
    M <- ncol(b_joint_std)
    J <- nrow(b_joint_std)
    b_type <- "std"
  }

  if(is.null(R_LD)){
    if(b_type == "non_std"){
      af <- check_af(af, J)
      sx <- sqrt(2*af*(1-af))
      b_joint_std <- b_joint*sx
    }
    if(!full_mat) return(colSums(b_joint_std^2))
    return(t(b_joint_std) %*% b_joint_std)
  }

  ld_mat <- check_R_LD(R_LD, "matrix")
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
  if(b_type == "non_std"){
    b_joint_std <- b_joint*sx
  }
  nb <- length(block_info$l)
  start_ix <- cumsum(c(1, block_info$l[-nb]))
  end_ix <- start_ix + block_info$l-1

  if(!full_mat){
    true_h2 <- lapply(seq(nb), function(i){
      # quad.diag(R, b_joint_std)
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
