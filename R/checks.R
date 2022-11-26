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

check_N <- function(N, n){
  if("matrix" %in% class(N) | "data.frame" %in% class(N)){
    if(ncol(N) ==1){
      N <- check_scalar_or_numeric(N, "N", n)
      Nc <- diag(n)
    }else{
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

check_R_LD <- function(R_LD){
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
  if(all(cl == "matrix")){
    R_LD <- lapply(R_LD, function(x){eigen(x)})
  }else if(any(cl == "matrix")){
    ii <- which(cl == "matrix")
    R_LD[ii] <- lapply(R_LD[ii], function(x){eigen(x)})
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


calc_ld_score <- function(R_LD){
  l2 <- lapply(R_LD, function(r){
    r2 <- with(r, tcrossprod(vectors, tcrossprod(vectors, diag(values)))^2)
    colSums(r2)
  }) %>% unlist()
  return(l2)
}
