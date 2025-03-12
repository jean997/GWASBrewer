check_scalar_or_numeric <- function(x, string, n){
  if(is.null(x)) return(x)
  if("matrix" %in% class(x) | "data.frame" %in% class(x)){
    if(ncol(x) > 1) stop(paste0(string, " must be a numeric vector or one column array."))
    x <- as.numeric(x[,1])
  }else if(length(x) == 1 & ("numeric" %in% class(x) | all(is.na(x)))){
    x <- rep(x, n)
  }else if(! "numeric" %in% class(x)){
    stop(paste0(string, " must be a numeric vector or one column array."))
  }
  if(length(x) != n) stop(paste0("Expected ", string, " to have length ", n, ", found ", length(x), "\n"))

  return(x)
}

check_matrix <- function(x, string, n, p){
  if(is.null(x)) return(x)
  if("data.frame" %in% class(x)){
    cat("Coercing ", string, " to matrix.\n")
    x <- as.matrix(x)
    colnames(x) <- NULL
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
  pi <- check_01(pi, "pi")
  return(pi)
}

check_N <- function(N, n, allow_mat = TRUE){
  if("sample_size" %in% class(N)){
    if(length(N$N) == n) return(N)
    stop(paste0("Expected information for ", n, " traits but found ", length(N$N), "\n"))
  }
  if(any(is.na(N))){
    stop("N cannot contain missing (NA) values. Use sample size 0 instead.\n")
  }
  if("matrix" %in% class(N) | "data.frame" %in% class(N)){
    if(ncol(N) == 1){
      type <- "vector"
    }else if("data.frame" %in% class(N) & "trait_1" %in% names(N)){
      type <- "df"
    }else{
      if(!allow_mat){
        stop("Matrix format for N is not allowed. Use scalar, vector, or data frame format.\n")
      }
      type <- "matrix"
    }
  }else{
    type <- "vector"
  }
  if(type == "df"){
    nn <- check_Ndf(N, n)
    return(nn)
  }
  if(type == "vector"){
    N <- check_scalar_or_numeric(N, "N", n)
    Nc <- diag(n)
    if(any(N == 0)){
      i <- which( N == 0)
      diag(Nc)[i] <- 0
    }
  }else if(type == "matrix"){
    N <- check_matrix(N, "N", n, n)
    max_val <- apply(N, 1, max)
    if(!all(max_val == diag(N)) | !all(N >=0) | !all.equal(N, t(N))){
      stop("N is not a valid sample size matrix.\nCheck that matrix is symmetric with only positive elements and that off-diagonal and diagonal entries are compatible.\n")
    }
    if(any(diag(N) == 0)){
      i <- which(diag(N) == 0)
      Nc <- matrix(0, nrow = n, ncol = n)
      Nci <- N[-i, -i, drop = FALSE] |> check_psd("N") |> stats::cov2cor()
      Nc[-i, -i] <- Nci
    }else{
      Nc <- check_psd(N, "N") |> stats::cov2cor()
    }
    N <- diag(N)
  }
  overlap <- !Matrix::isDiagonal(Nc)
  ret <- list("Nc"= Nc, "N" = N, "overlap" = overlap) |> structure(class = c("sample_size", "list"))
  return(ret)
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
  if(!all(names(N) %in% c(paste0("trait_", 1:M), "N"))){
    unused_cols <- setdiff(names(Ndf), c(paste0("trait_", 1:M), "N"))
    warning(paste0("Sample size dataframe contains unused columns: ", paste0(unused_cols, collapse = ","), "\n"))
  }
  N <- dplyr::select(N, all_of(c(paste0("trait_", 1:M), "N"))) %>%
        mutate_at(paste0("trait_", 1:M), as.logical)
  nr <- nrow(N)
  nd <- select(N, -N) %>% distinct() %>% nrow()
  if(nr != nd){
    warning(paste0("Sample size dataframe contains non-unique study combinations. Duplicated rows will be collapsed\n"))
    N <- N %>% group_by_at(paste0("trait_", 1:M)) %>% summarize(N = sum(N)) %>% ungroup()
  }
  trait_present <- N %>% select(-N) %>% apply(2, max)
  if(!all(trait_present)){
    warning(paste0("Not all traits have sample size given in data frame N. Missing ", paste0(which(!trait_present), collapse = ","), "\n"))
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
  ct <- diag(Nmat)
  if(any(ct == 0)){
    i <- which(ct == 0)
    Nc <- matrix(0, nrow = M, ncol = M)
    Nc[-i, -i] <- stats::cov2cor(Nmat[-i, -i])
  }else{
    Nc <- stats::cov2cor(Nmat)
  }
  overlap <- !Matrix::isDiagonal(Nc)
  ret <- list("Ndf" = N, "Nc"= Nc, "N" = ct, "overlap" = overlap) |> structure(class = c("sample_size", "sample_size_df", "list"))
  return(ret)
}

subset_N_nonzero <- function(N){
  stopifnot("sample_size" %in% class(N))
  i <- which(N$N > 0)
  newN <- N
  newN$N <- N$N[i]
  newN$Nc <- N$Nc[i,i]
  newN$zero_ix <- setdiff(1:length(N$N), i)
  newN$nonzero_ix <- i
  return(newN)
}



check_psd <- function(M, string, tol = 1e-8){
  #if(!Matrix::isSymmetric(M)){
  #  stop(paste0(string, " must be symmetric.\n"))
  #}
  ## all.equal uses sqrt(.Machine$double.eps) as tolerance by default
  if(!all.equal(M, t(M))){
    stop(paste0(string, " must be symmetric.\n"))
  }
  eMvals <- fast_eigen_vals(M)
  if(!all(eMvals >= -1*tol)){
    stop(paste0(string, " is not positive semi-definite.\n"))
  }
  return(M)
}


check_G <- function(G, h2){
  if(! "matrix" %in% class(G)){
    if(!(class(G) == "numeric" | class(G) == "integer" )){
      stop(paste0("G should have class matrix, numeric, or integer, found ", class(G), "\n"))
    }
    if(!(G >= 0 & G == round(G))){
      stop("If G is not a matrix, it should be a non-negative integer. Found ", G, "\n")
    }
    if(G == 0 | G == 1){
      n <- 1
    }else{
      n <- G
    }
    G <- matrix(0, nrow = n, ncol = n)
  }else{
    n <- nrow(G)
  }


  G <- check_matrix(G, "G", n, n)
  h2 <- check_scalar_or_numeric(h2, "h2", n)
  h2 <- check_01(h2, "h2")
  if(!all(diag(G) ==0)){
    stop("G must have 0s on the diagonal.")
  }
  G_tot <- tryCatch(direct_to_total(G), error = function(e){
    stop("Failed to compute total effects from direct. Check that supplied G corresponds to a valid DAG.\n")
  })
  if(!all(diag(G_tot) == 0)){
    stop("Supplied G does not correspond to a valid DAG.\n")
  }
  direct_h2_1 <- solve(diag(n) + t(G_tot)^2) %*% rep(1, n)
  if(any(direct_h2_1 < 0)){
    stop("Supplied G does not correspond to a valid DAG.\n")
  }
  direct_h2 <- solve(diag(n) + t(G_tot)^2) %*% h2
  direct_e2 <- solve(diag(n) + t(G_tot)^2) %*% (1-h2)
  if(any(direct_h2 < 0)){
    stop(paste0("Supplied G is incompatible with supplied h2. You could try increasing the heritability of traits ",
                paste0(which(direct_h2 < 0), collapse = ","), ".\n"))
  }
  return(list(G_dir = G, G_tot = G_tot, dir_h2 = as.vector(direct_h2), h2 = h2, M = n, dir_e2 = as.vector(direct_e2)))
}

check_01 <- function(x, name){
  if(any(x <0) | any(x > 1)) stop(paste0("Expected ", name, " to be between 0 and 1. Check elements",
                                    paste0(which(x < 0 | x > 1), collapse = ",")))
  return(x)
}

direct_to_total <- function(G_dir){
  n <- nrow(G_dir)
  G_total <- G_dir %*% solve(diag(n) - G_dir)
  if(!all(diag(G_total) == 0)){
    stop("Failed to compute total effects from direct. Check that supplied G corresponds to a valid DAG.\n")
  }
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
      R_LD <- lapply(R_LD, function(x){fast_eigen(x)})
    }else if(any(cl == "matrix")){
      ii <- which(cl == "matrix")
      R_LD[ii] <- lapply(R_LD[ii], function(x){fast_eigen(x)})
    }
  }else if(return=="matrix"){
    if(all(cl == "eigen")){
      R_LD <- lapply(R_LD, function(x){
        tcrossprod(x$vectors, (x$vectors * rep(x$values, each = nrow(x$vectors))))
        #tcrossprod(x$vectors, tcrossprod(x$vectors, diag(x$values)))
      })
    }else if(any(cl == "eigen")){
      ii <- which(cl == "eigen")
      R_LD[ii] <- lapply(R_LD[ii], function(x)function(x){
        tcrossprod(x$vectors, (x$vectors * rep(x$values, each = nrow(x$vectors))))
        #tcrossprod(x$vectors, tcrossprod(x$vectors, diag(x$values)))
      })
    }else{
      R_LD <- lapply(R_LD, as.matrix)
    }
  }else if(return == "sqrt"){
    if(all(cl == "eigen")){
      R_LD <- lapply(R_LD, function(x){
        x$vectors*rep(sqrt(x$values), each = nrow(x$vectors))
        #tcrossprod(x$vectors, diag(sqrt(x$values)))
      })
    }else if(all(cl == "matrix")){
      R_LD <- lapply(R_LD, function(m){
        x <- fast_eigen(m)
        x$vectors*rep(sqrt(x$values), each = nrow(x$vectors))
        #tcrossprod(x$vectors, diag(sqrt(x$values)))
      })
    }else{
      ie <- which(cl == "eigen")
      R_LD[ie] <- lapply(R_LD[ie], function(x){
        x$vectors*rep(sqrt(x$values), each = nrow(x$vectors))
        #tcrossprod(x$vectors, diag(sqrt(x$values)))
      })
      im <- which(cl == "matrix")
      R_LD[im] <- lapply(R_LD[im], function(m){
        x <- fast_eigen(m)
        x$vectors*rep(sqrt(x$values), each = nrow(x$vectors))
        #tcrossprod(x$vectors, diag(sqrt(x$values)))
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



check_snpinfo <- function(snp_info, J){
  if(!"data.frame" %in% class(snp_info)){
    stop("snp_info should have class containing data.frame.\n")
  }
  if(nrow(snp_info) != J){
    stop(paste0("Expected snp_info to have ", J, " rows but found ", nrow(snp_info), "."))
  }
  if(any(c("SNP", "AF") %in% names(snp_info))){
    warning("Provided snp_info contains columns SNP or AF. These will be over-written.")
  }
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
  af <- check_01(af, "af")
  return(af)
}

check_effect_function_list <- function(snp_effect_function, M, snp_info = NULL){
  if(!class(snp_effect_function) == "list"){
    f <- check_snp_effect_function(snp_effect_function, snp_info)
    fl <- list()
    for(i in 1:M) fl[[i]] <- f
    return(fl)
  }else{
    if(length(snp_effect_function) != M){
      stop(paste0("Expected snp_effect_function to be a list of length ", M, ", a function, or 'normal'"))
    }
    fl <- list()
    for(i in 1:M) fl[[i]] <- check_snp_effect_function(snp_effect_function[[i]], snp_info)
    return(fl)
  }
}

check_snp_effect_function <- function(snp_effect_function, snp_info = NULL){

  if(!(class(snp_effect_function) == "character" | class(snp_effect_function) == "function") ){
    stop("snp_effect_function should either be 'normal', or a function.")
  }
  if(class(snp_effect_function) == "function"){
    if(!is.null(snp_info)){
      if(nrow(snp_info) < 1000 ){
        nr <- ceiling(1000/nrow(snp_info))
        test_snp_info <- do.call("rbind", replicate(nr, snp_info, simplify = FALSE))[1:1000,]
        test_snp_info$SNP <- 1:1000
      }else{
        test_snp_info <- snp_info[1:1000,]
      }
    }else{
      test_snp_info <- data.frame(SNP = 1:1000, AF = stats::rbeta(n = 1000, 1, 5))
    }
    x <- tryCatch(snp_effect_function(n = 1000, sd = 1, snp_info = test_snp_info), error = function(e){
      stop(paste0("Failed to run snp_effect_function with error: ", e) )
    })
    test_stat <- try(stats::t.test(x^2 - (1/1000)), silent = TRUE)
    if(!inherits(test_stat, "try-error")){
      if(test_stat$p.value < 0.01){
        warning(paste0("Your snp_effect_function may not generate effects with the correct variance. This
                will lead to realized heritabilities that are very different from your desired heritability.
                This warning was generated in a test run. I generated 1000 effect sizes from your function
                and expected to find the total sum of squares approximately equal to 1. The total sum of
                squares I observed was ", sum(x^2), ". This warning will occur occaisionally for correct
                functions due to randomness. If you see this warning, verify that your function is correct.
                Check the vignette or ?sim_mv for more details."))
      }
    }
    return(snp_effect_function)
  }else if(snp_effect_function == "normal"){
    f <- function(n, sd, ...){
      stats::rnorm(n =n, mean = 0, sd = sd/sqrt(n))
    }
    return(f)
  }
}


# calc_ld_score <- function(R_LD){
#   l2 <- lapply(R_LD, function(r){
#     r2 <- with(r, tcrossprod(vectors, tcrossprod(vectors, diag(values)))^2)
#     colSums(r2)
#   }) %>% unlist()
#   return(l2)
# }

#https://stackoverflow.com/questions/27870542/r-check-if-function-in-a-package-was-called-from-the-package-fun-or-externally
called_intern <- function() {
  te <- topenv(parent.frame(1))
  if(isNamespace(te) && getNamespaceName(te) == "GWASBrewer") {
    return(1) # called internally
  }
  return(0)
}
