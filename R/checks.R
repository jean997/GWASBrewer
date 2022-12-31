check_scalar_or_numeric <- function(x, string, n){
  if(is.null(x)) return(x)
  if("matrix" %in% class(x) | "data.frame" %in% class(x)){
    if(ncol(x) > 1) stop(paste0(string, " must be a numeric vector or one column array."))
    x <- as.numeric(x[,1])
  }else if(! "numeric" %in% class(x)){
    stop(paste0(string, " must be a numeric vector or one column array."))
  }else if(length(x) == 1){
    x <- rep(x, n)
  }
  if(length(x) != n) stop(paste0("Expected ", string, " to have length ", n, ", found ", length(x), "\n"))

  return(x)
}

check_matrix <- function(x, string, n, p){
  if(is.null(x)) return(x)
  if("data.frame" %in% class(x)){
    cat("Coercing ", string, " to matrix.\n")
    x <- as.matrix(x)
  }else if("numeric" %in% class(x)){
    x <- as.matrix(x, ncol = 1)
  }else if(!"matrix" %in% class(x)){
    stop(paste0(string, " must be a numeric vector, matrix, or data.frame."))
  }
  if(!missing(n)){
    if(nrow(x) != n) stop(paste0("Expected ", string, " to have ", n, " rows, found ", nrow(x), "\n"))
  }
  if(!missing(p)){
    if(ncol(x) != p) stop(paste0("Expected ", string, " to have ", p, " columns, found ", ncol(x), "\n"))
  }
  return(x)
}


check_pi <- function(pi, n, p){
  if("matrix" %in% class(pi)){
    pi <- check_matrix(pi, "pi", n, p)
  }else{
    pi <- check_scalar_or_numeric(pi, "pi", p)
  }
  pi <- check_01(pi)
  return(pi)
}

check_N <- function(N, n, allow_mat = TRUE){
  if("matrix" %in% class(N) | "data.frame" %in% class(N)){
    if(ncol(N) ==1){
      N <- check_scalar_or_numeric(N, "N", n)
      Nc <- diag(n)
    }else{
      if("data.frame" %in% class(N) & "trait_1" %in% names(N)){
        nn <- check_Ndf(N, n)
        return(nn)
      }
      if(!allow_mat){
        stop("Matrix format for N is not allowed. Use scalar, vector, or data frame format.\n")
      }
      N <- check_matrix(N, "N", n, n)
      N <- check_psd(N, "N")
      Nc <- cov2cor(N)
      N <- diag(N)
    }
  }else{
    N <- check_scalar_or_numeric(N, "N", n)
    Nc <- diag(n)
  }
  return(list("Nc"= Nc, "N" = N))
}

make_Ndf_indep <- function(N){
  M <- length(N)
  df <- data.frame(diag(M))
  names(df) <- paste0("trait_", 1:M)
  df$N <- N
  df <- mutate_at(df,paste0("trait_", 1:M), as.logical)
  return(df)
}

check_Ndf <- function(N, M){
  if(!"data.frame" %in% class(N)){
    stop("class(N) does not include data.frame\n")
  }
  if(missing(M)){
    M <- ncol(N)-1
  }
  if(!"N" %in% names(N)){
    stop("Did not find column N in data frame N")
  }
  if(!all(paste0("trait_", 1:M) %in% names(N))){
    stop(paste0("Did not find all of trait_1 ...trait_", M, " in data frame N\n"))
  }
  N <- dplyr::select(N, all_of(c(paste0("trait_", 1:M), "N"))) %>%
        mutate_at(paste0("trait_", 1:M), as.logical)
  nr <- nrow(N)
  nd <- select(N, -N) %>% distinct() %>% nrow()
  if(nr != nd){
    warning(paste0("Data frame N contains non-unique study combinations. Duplicated rows will be collapsed\n"))
    N <- N %>% group_by_at(paste0("trait_", 1:M)) %>% summarize(N = sum(N)) %>% ungroup()
  }
  trait_present <- N %>% select(-N) %>% apply(2, max)
  if(!all(trait_present)){
    stop(paste0("Not all traits have sample size given in data frame N. Missing ", paste0(which(!trait_present), collapse = ","), "\n"))
  }


  Nmat_long <- expand.grid(1:M, 1:M) %>% filter(Var1 <= Var2)
  #Nmat_long$N <- NA
  Nmat_long$N <- sapply(seq(nrow(Nmat_long)), function(i){
    v1 <- paste0("trait_", Nmat_long$Var1[i])
    v2 <- paste0("trait_", Nmat_long$Var2[i])
    filter_at(N, vars(all_of(c(v1, v2))), all_vars(. == TRUE)) %>% with(., sum(N))
  })
  Nmat <- reshape2::acast(Nmat_long, Var1 ~ Var2, value.var = "N", drop = FALSE, fill = 0)
  Nmat2 <-reshape2::acast(Nmat_long, Var2 ~ Var1, value.var = "N", drop = FALSE, fill = 0)
  diag(Nmat2) <- 0
  Nmat <- Nmat + Nmat2
  Nc <- cov2cor(Nmat)
  return(list(Ndf = N, N = diag(Nmat), Nc = Nc))
}

check_psd <- function(M, string){
  if(!Matrix::isSymmetric(M)){
    stop(paste0(string, " must be symmetric.\n"))
  }
  eM <- eigen(M)
  if(!all(eM$vaues >= 0)){
    stop(paste0(string, " is not positive semi-definite.\n"))
  }
  return(M)
}

check_G <- function(G, h2, n){
  G <- check_matrix(G, "G", n, n)
  h2 <- check_scalar_or_numeric(h2, "h2", n)
  if(!all(diag(G) ==0)){
    stop("G must have 0s on the diagonal.")
  }
  G_tot <- try(direct_to_total(G), silent = TRUE)
  if("try-error" %in% class(G_tot)){
    stop("Failed to compute total effects from direct. Check that supplied G corresponds to a valid DAG.\n")
  }
  if(!all(diag(G_tot) == 0)){
    stop("Supplied G does not correspond to a valid DAG.\n")
  }
  direct_h2_1 <- solve(diag(n) + t(G_tot)^2) %*% rep(1, n)
  if(any(direct_h2_1 < 0)){
    stop("Supplied G does not correspond to a valid DAG.\n")
  }
  direct_h2 <- solve(diag(n) + t(G_tot)^2) %*% h2
  if(any(direct_h2 < 0)){
    stop(paste0("Supplied G is incompatible with supplied h2. You could try increasing the heritability of traits ",
                paste0(which(direct_h2 < 0), collapse = ","), ".\n"))
  }
  return(list(G_dir = G, G_tot = G_tot, dir_h2 = as.vector(direct_h2)))
}

check_01 <- function(x, name){
  if(any(x <0) | any(x > 1)) stop(paste0("Expected ", name, " to be between 0 and 1. Check elements",
                                    paste0(which(x < 0 | x > 1), collapse = ",")))
  return(x)
}

direct_to_total <- function(G_dir){
  n <- nrow(G_dir)
  G_total <- G_dir %*% solve(diag(n) - G_dir)
  return(G_total)
}

check_R_LD <- function(R_LD, return = c("eigen", "matrix", "sqrt", "l")){
  return <- match.arg(return, return)
  if(class(R_LD) != "list"){
    stop(paste0("R_LD should be of class list, found ", class(R_LD), "\n"))
  }
  cl <- sapply(R_LD, function(x){
    case_when("matrix" %in% class(x) ~ "matrix",
              class(x) == "dsCMatrix" ~ "matrix",
              class(x) == "eigen" ~ "eigen",
              TRUE ~ "not_allowed")
  })
  if(any(cl == "not_allowed")){
    stop("R_LD should be a list with elements of class matrix, dsCMatrix, or eigen.")
  }
  if(return=="eigen"){
    if(all(cl == "matrix")){
      R_LD <- lapply(R_LD, function(x){eigen(x)})
    }else if(any(cl == "matrix")){
      ii <- which(cl == "matrix")
      R_LD[ii] <- lapply(R_LD[ii], function(x){eigen(x)})
    }
  }else if(return=="matrix"){
    if(all(cl == "eigen")){
      R_LD <- lapply(R_LD, function(x){
        tcrossprod(x$vectors, tcrossprod(x$vectors, diag(x$values)))
      })
    }else if(any(cl == "eigen")){
      ii <- which(cl == "eigen")
      R_LD[ii] <- lapply(R_LD[ii], function(x)function(x){
        tcrossprod(x$vectors, tcrossprod(x$vectors, diag(x$values)))
      })
    }else{
      R_LD <- lapply(R_LD, as.matrix)
    }
  }else if(return == "sqrt"){
    if(all(cl == "eigen")){
      R_LD <- lapply(R_LD, function(x){
        tcrossprod(x$vectors, diag(sqrt(x$values)))
      })
    }else if(all(cl == "matrix")){
      R_LD <- lapply(R_LD, function(m){
        x <- eigen(m)
        tcrossprod(x$vectors, diag(sqrt(x$values)))
      })
    }else{
      ie <- which(cl == "eigen")
      R_LD[ie] <- lapply(R_LD[ie], function(x){
        tcrossprod(x$vectors, diag(sqrt(x$values)))
      })
      im <- which(cl == "matrix")
      R_LD[im] <- lapply(R_LD[im], function(m){
        x <- eigen(m)
        tcrossprod(x$vectors, diag(sqrt(x$values)))
      })
    }
  }else if(return == "l"){
    if(all(cl == "eigen")){
      l <- sapply(R_LD, function(x){length(x$values)})
    }else if(all(cl == "matrix")){
      l <- sapply(R_LD, nrow)
    }else{
      l <- rep(NA, length(R_LD))
      ie <- which(cl == "eigen")
      l[ie] <- sapply(R_LD[ie], function(x){length(x$values)})
      im <- which(cl == "matrix")
      l[im] <- sapply(R_LD[im], function(m){nrow(m)})
    }
    return(l)
  }
  return(R_LD)
}


# l is list of block lengths
check_snpinfo <- function(snp_info, l){
  if(nrow(snp_info) != sum(l)){
    stop(paste0("R_LD contains information for ", sum(l), " variants but snp_info contains information for ", nrow(snp_info), " variants."))
  }
  if(!all(c("SNP", "AF") %in% names(snp_info))){
    stop("snp_info must contain columns AF and SNP")
  }
  snp_info$block <- rep(seq(length(l)), l)
  snp_info$ix_in_block <- sapply(l, function(ll){seq(ll)}) %>% unlist()
  return(snp_info)
}

check_af <- function(af, n, function_ok = TRUE){
  if(is.null(af)){
    return(NULL)
  }else if(class(af) == "function"){
    if(!function_ok) stop("af cannot be a function.\n")
    myaf <- af(n)
    af <- myaf
  }
  af <- check_scalar_or_numeric(af, "af", n)
  af <- check_01(af)
  return(af)
}



check_snp_effect_function <- function(snp_effect_function){

  if(!(class(snp_effect_function) == "character" | class(snp_effect_function) == "function")){
    stop("snp_effect_function should either be 'normal' or a function.")
  }
  if(class(snp_effect_function) == "function"){

    x <- tryCatch(snp_effect_function(n = 3, sd = 1), error = function(e){
      stop(paste0("Failed to run snp_effect_function with error: ", e) )
    })
    x <- tryCatch(snp_effect_function(n = 3, sd = 1,af= rep(0.5, 3)), error = function(e){
      stop(paste0("Failed to run snp_effect_function with error: ", e) )
    })
    return(snp_effect_function)
  }else if(snp_effect_function == "normal"){
    f <- function(n, sd, ...){
      rnorm(n =n, mean = 0, sd = sd)
    }
    return(f)
  }else{
    stop("snp_effect_function should either be 'normal' or a function.")
  }
}


calc_ld_score <- function(R_LD){
  l2 <- lapply(R_LD, function(r){
    r2 <- with(r, tcrossprod(vectors, tcrossprod(vectors, diag(values)))^2)
    colSums(r2)
  }) %>% unlist()
  return(l2)
}
