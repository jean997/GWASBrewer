
#'@title Sample individual level data with joint effects matching a sim_mv object
#'@param N Sample size, scalar, vector, or special sample size format data frame, see details.
#'@param dat An object of class \code{sim_mv} (produced by \code{sim_mv}). If `dat` is omitted,
#'the function will generate a matrix of genotypes only. If `dat` is provided, phenotypes for the traits
#'in `dat` will also be included.
#'@param genos Optional matrix of pre-generated genotypes. If \code{genos} is supplied, \code{resample_inddata} will only generate phenotypes.
#'@param J Optional number of variants. \code{J} is only required if \code{dat} is missing.
#'@param R_LD LD pattern (optional). See \code{?sim_mv} for more details.
#'@param af Allele frequencies. \code{af} is required unless unless \code{genos} is supplied.
#'@param new_env_var Optional. The environmental variance in the new population.
#'If missing the function will assume the environmental variance is the same as in the old population.
#'@param new_R_E Optional, specify environmental correlation in the new population.
#'If missing, the function will assume the environmental correlation is the same as in the original data.
#'@param new_R_obs Optional, specify overall trait correlation in the new population. Specify at most one of \code{new_R_E} or \code{new_R_obs}.
#'If missing, the function will assume the environmental correlation is the same as in the original data.
#'@param calc_sumstats If \code{TRUE}, associations between genotypes and phenotypes will be calculated and returned.
#'@details This function can be used to generate individual level genotype and phenotype data. It can be used in three modes:
#'
#'To generate genotype data only: No \code{sim_mv} object needs to be included. Supply only \code{N} as a single integer for the
#'number of individuals, \code{J} for the number of variants, \code{af}, and \code{R_LD} if desired. All other
#'parameters are not relevant if there is no phenotype, so if they are supplied, you will get an error. The returned object will include a
#'\code{N x J} matrix of genotypes and a vector of allele frequencies.
#'
#'To generate both genotype and phenotype data: Supply \code{dat} (a \code{sim_mv} object) and leave \code{genos} missing. \code{N} and \code{af}
#'are required and all other options are optional.
#'
#'To generate phenotype data only: Supply \code{dat} (a \code{sim_mv} object) and provide a matrix of genotypes to the \code{genos} argument. The number of
#'rows in \code{genos} must be equal to the total number of individuals implied by \code{N}.
#' So for example, if there are two traits with 10 samples each and no overlap, \code{genos}
#'should have 20 rows. The \code{R_LD} and \code{af} arguments should contain the population
#'LD and allele frequencies used to produce the genotypes. These are used to compute the genetic variance-covariance matrix.
#'\code{N} and \code{af} are required and all other options are optional.
#'
#'@examples
#' # Use gen_gwas_from_b to generate individual level data with given effect size.
#' Ndf <- data.frame(trait_1 = 1, trait_2 = 1, N = 10000)
#' G <- matrix(0, nrow = 2, ncol = 2)
#' R_E <- matrix(c(1, 0.8, 0.8, 1), nrow = 2, ncol = 2)
#' # original data
#' dat <- sim_mv(N = Ndf, J = 2000, h2 = c(0.4, 0.3), pi = 100/2000,
#'                G = G, R_E = R_E, af = function(n){rbeta(n, 1, 5)})
#'# Now generate GWAS data
#'gw_dat <- resample_inddata(dat, N = Ndf, calc_sumstats = TRUE)
#'@export
resample_inddata <- function(N,
                             dat = NULL,
                             genos = NULL,
                             J = NULL,
                             R_LD = NULL,
                             af = NULL,
                             sim_func = gen_genos_mvn,
                             new_env_var = NULL,
                             new_R_E = NULL,
                             new_R_obs = NULL,
                             calc_sumstats = FALSE){

  if(is.null(af)){
    stop("af is required to generate individual level data.\n")
  }

  # Option 1: generate genotypes and phenotypes provide dat, omit genos and J
  # Option 2: generate genotype data only. dat is missing, provide N, J, af, and possibly R_LD
  # Option 3: Generate phenotypes for existing genotypes. Provide genos, and dat. omit J

  if(is.null(dat)){
    message("Generating genotype matrix only.")
    if(! "numeric" %in% class(N)| "integer" %in% class(N) | length(N) > 1 | ! N == round(N) | N <= 0){
      stop("If dat is missing, N should be a positive integer.")
    }
    M <- 1
    if(!is.null(new_env_var) | !is.null(new_R_E) | !is.null(new_R_obs) | calc_sumstats){
      stop("If dat is omitted, new_env_var, new_R_E, new_R_obs should be omitted and calc_sumstats should be FALSE. See help page for more details.\n")
    }
    if(!is.null(genos)){
      stop("genos were provided so I have nothing to do.\n")
    }
    if(is.null(J)){
      stop("Provide J if dat is missing.\n")
    }
  }else{
    if(!"sim_mv" %in% class(dat)) stop("dat should have class sim_mv.")
    if(!is.null(J)){
      stop("If dat is provided, don't provide J.")
    }
    M <- ncol(dat$beta_hat)
    J <- nrow(dat$beta_hat)
    if(!is.null(genos)){
      #if(!is.null(R_LD) | !is.null(af)){
      #  warning("Genotypes are provided so supplied R_LD and af are ignored.")
      #}
      message("Generating phenotypes only.")
    }else{
      message("Generating both genotypes and phenotypes.")
    }
  }



  ## check sample size and generate genotype data
  nn <- check_N(N, M, allow_mat = FALSE)
  if(is.null(nn$Ndf)){
    nn$Ndf <- make_Ndf_indep(nn$N)
  }
  ntotal <- sum(nn$Ndf$N)

  ## simulate genotypes
  if(is.null(genos)){
    if(is.null(R_LD)){
      af <- check_af(af, J)
      X <- sim_func(ntotal, J, NULL, af)
    }else{
      l <- check_R_LD(R_LD, "l")
      af <- check_af(af, sum(l), function_ok = FALSE)
      X <- sim_func(ntotal, J, R_LD, af)
    }
    my_af <- check_af(X$af, J, function_ok = FALSE)
    X <- check_matrix(X$X, "genos", ntotal, J)
    if(is.null(dat)){
      R <- list(X = X,af = my_af)
      return(R)
    }
  }else{
    X <- check_matrix(genos, "genos", ntotal, J)
    if(! all(X %in% 0:2)){
      stop("All elements of genos should be 0, 1, or 2.")
    }
  }

  message(paste0("SNP effects provided for ", J, " SNPs and ", M, " traits."))

  if( dat$geno_scale == "sd"){
    message("Original data have effects on the per-genotype sd scale. I will assume that per-genotype sd effects are the same in the new and old populations.")
    # Compute non-standardized effects
    nr <- ceiling(J/length(af))
    dat <- rescale_sumstats(dat = dat, output_geno_scale = "allele", af = rep(af, nr)[1:J],
                            output_pheno_sd = dat$pheno_sd)

  }
  # compute variance explained by genetics in new population
  # not heritability despite using compute_h2 function
  Sigma_G <- compute_h2(b_joint = dat$beta_joint,
                        geno_scale = dat$geno_scale,
                        pheno_sd = 1,
                        R_LD = R_LD,
                        af = af,
                        full_mat = TRUE)

  v_G <- diag(Sigma_G)
  if(!all(Sigma_G == dat$Sigma_G)){
    message("Genetic variance in the new population differs from the genetic variance in the old population.")
  }

  ## Environmental variance
  if(is.null(new_env_var)){
    message("I will assume that the environmental variance is the same in the old and new population.")
    v_E <- diag(dat$Sigma_E)
  }else{
    v_E <- check_scalar_or_numeric(new_env_var, "new_env_var", M)
  }
  pheno_sd <- sqrt(v_G + v_E)
  h2 <- v_G/(v_G + v_E)
  if(is.null(new_R_obs) & is.null(new_R_E)){
    message("I will assume that environmental correlation is the same in the old and new population. Note that this could result in different overall trait correlations.")
    R_E <- cov2cor(dat$Sigma_E)
    Sigma_E <- diag(sqrt(v_E),nrow = M) %*% R_E %*% diag(sqrt(v_E), nrow = M)
    trait_corr <- cov2cor(Sigma_G + Sigma_E)
  }else if(!is.null(new_R_obs)){
    new_R_obs <- check_matrix(new_R_obs, "new_R_obs", M, M)
    new_R_obs <- check_psd(new_R_obs, "new_R_obs")
    if(!all(diag(new_R_obs) == 1)) stop("new_R_obs should be a correlation matrix. Found diagonal entries not equal to 1.")
    trait_corr <- new_R_obs
    Sigma_tot <- diag(pheno_sd, nrow  = M) %*% new_R_obs %*% diag(pheno_sd, nrow = M)
    Sigma_E <- Sigma_tot - Sigma_G
    Sigma_E <- tryCatch(check_psd(Sigma_E, "Sigma_E"), error = function(e){
      stop("new_R_obs is incompatible with trait relationships and heritability.")
    })
  }else if(!is.null(new_R_E)){
    new_R_E <- check_matrix(new_R_E, "new_R_E", M, M)
    new_R_E <- check_psd(new_R_E, "new_R_E")
    if(!all(diag(new_R_E) == 1)) stop("new_R_E should be a correlation matrix. Found diagonal entries not equal to 1.")
    Sigma_E <- diag(sqrt(v_E),nrow = M) %*% new_R_E %*% diag(sqrt(v_E), nrow = M)
    trait_corr <- cov2cor(Sigma_G + Sigma_E)
  }
  if(!all(pheno_sd == dat$pheno_sd)){
    message("Note that the phenotype in the new population has a different variance from the phenotype in the old population.")
    message("I will keep the phenotype on the same scale as the original data, so effect sizes in the old and new object are comparable.")
  }

  eV <- eigen(Sigma_E)

  study_id <- rep(1:nrow(nn$Ndf), nn$Ndf$N)
  study_ix <- lapply(1:M, function(m){
    ix <- which(nn$Ndf[paste0("trait_", m)] == TRUE)
    which(study_id %in% ix)
  })

  # environmental component of phenotype
  y_E <- matrix(rnorm(n = M*ntotal), nrow = ntotal)
  y_E <- y_E %*% diag(sqrt(eV$values), M) %*% t(eV$vectors)

  # YG <- purrr::map_dfc(1:M, function(m){
  #   b <- dat$beta_joint[,m]
  #   y_gen <- t(t(X[study_ix[[m]],])*b) %>% rowSums()
  #   my_y <- data.frame(y = rep(NA, ntotal))
  #   names(my_y)[1] <- paste0("y_", m)
  #   my_y[study_ix[[m]],1] <- y_gen
  #   return(my_y)
  # })
  # YE <- purrr::map_dfc(1:M, function(m){
  #   my_y <- data.frame(y = rep(NA, ntotal))
  #   names(my_y)[1] <- paste0("y_", m)
  #   my_y[study_ix[[m]],1] <- y_E[study_ix[[m]], m]
  #   return(my_y)
  # })
  Y <- purrr::map_dfc(1:M, function(m){
    b <- dat$beta_joint[,m]
    y_gen <- t(t(X[study_ix[[m]],])*b) %>% rowSums()
    y <- y_gen + y_E[study_ix[[m]], m]
    my_y <- data.frame(y = rep(NA, ntotal))
    names(my_y)[1] <- paste0("y_", m)
    my_y[study_ix[[m]],1] <- y
    return(my_y)
  })


  R <- list(X = X,
            Y= Y,
            # YG = YG,
            # YE = YE,
            af = af,
            Sigma_G = Sigma_G,
            Sigma_E = Sigma_E,
            pheno_sd = pheno_sd,
            h2 = h2,
            trait_corr = trait_corr,
            beta_joint = dat$beta_joint)
  if(!is.null(genos)) R$X <- NULL
  if(!calc_sumstats) return(R)

  ## Calculate GWAS estimates fast
  sumstats <- fast_lm(X, Y, check = FALSE)
  R$beta_hat <- dplyr::select(sumstats, all_of(paste0("bhat_", 1:M)) ) %>% as.matrix()
  R$s_estimate <- dplyr::select(sumstats, all_of(paste0("s_", 1:M)) ) %>% as.matrix()
  return(R)

}

#'@export
fast_lm <- function(X, Y, check = TRUE){
  if(check){
    X <- check_matrix(X, "X")
    N <- nrow(X)
    J <- ncol(X)
    Y <- check_matrix(Y, "Y", N)
    M <- ncol(Y)
  }else{
    N <- nrow(X)
    J <- ncol(X)
    M <- ncol(Y)
  }
  sumstats <- purrr::map_dfc(1:M, function(m){
    x_ix <- which(!is.na(Y[,m]))
    yc <- Y[x_ix,m] - mean(Y[x_ix,m])
    xc <- scale(X[x_ix, ], center = TRUE, scale = FALSE)
    xty <- colSums(xc*yc)
    xtx <- colSums(xc^2)
    bhat <- xty/xtx
    yhat <- t(t(xc)*bhat)
    r2 <- colSums((yhat - yc)^2)
    s2 <- r2/(xtx*(length(yc)-2))
    df <- data.frame(bhat = bhat, s = sqrt(s2))
    names(df) <- paste0(c("bhat_", "s_"), m)
    return(df)
  })
  return(sumstats)
}


#'@export
gen_genos_mvn <- function(n, J, R_LD, af){

  if(is.null(R_LD)){
    af <- check_af(af, J)
    X <- replicate(n = n, rbinom(n = J, size = 2, prob = af)) %>% t()
    if(J == 1) X <- t(X)
    return(list(X = X, af = af))
  }
  # Checks happen in R_LD_to_haplodat
  l <- check_R_LD(R_LD, "l")
  #af <- check_af(af, sum(l), function_ok = FALSE)


  hdat <- R_LD_to_haplodat(R_LD, af = af)
  block_info <- assign_ld_blocks(l, J)
  nb <- length(l)
  if(!is.null(block_info$last_block_info)){
    b <- block_info$last_block_info[1]
    x <- block_info$last_block_info[2]
    last_block <- hdat[[b]]$cor[seq(x), seq(x)]
    last_af<- hdat[[b]]$freqs[seq(x)]
    hdat_last <- R_LD_to_haplodat(R_LD = list(last_block), last_af)
    hdat[[nb + 1]] <- hdat_last[[1]]
    l <- c(l, x)
    nb <- nb + 1
  }

  X <- purrr::map_dfr(block_info$block_index, function(bi){
    xx <- hapsim_simple(n = 2*n, hap = hdat[[bi]])
    a <- xx[seq(1, 2*n, by = 2), ] + xx[seq(2, 2*n, by = 2), ]
    a <- data.frame(t(a))
    names(a) <- paste0("S", 1:n)
    return(a)
  })

  return(list(X = t(X), af = af[block_info$index]))


}
