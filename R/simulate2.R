#'@title Simulate multivariate GWAS data
#'@param N Either a single number representing GWAS sample size for all studies or a
#'vector of length equal to the number of studies that will be generated.
#'@param J Number of variants to simulate
#'@param h2 A vector or single number giving the heritability of each trait.
#'@param pi A vector or single number giving the expected proportion of direct effect SNPs for each trait.
#'@param G If using 'general' mode (see details), G is a matrix of direct effects. Rows correspond to the 'from' trait
#'and columns correspond to the 'to' variable, so G[1,2] is the direct effect of variable 1 on variable 2.
#'@param taux_xz,tau_yz Used in 'xyz' mode (see details below)
#'@param dir_xz,dir_yz Used in 'xyz' mode (see details below)
#'@param gamma Used in 'xyz' mode (see details below)
#'@param R_E Environmental correlation
#'@param overlap_prop Proportion of GWAS samples overlapping between studies. Scalar.
#'@param R_LD List of eigen decompositions of LD correlation matrices, may be missing.
#'@param snp_info If R_LD is provided, provide a data frame with columns "SNP" and "AF"
#'
#'@return A list with the following elements:
#'
#'Simulated effect estimates and standard errors are contained in matrices
#'
#' + `beta_hat` Effect estimates for each trait
#'
#' + `se_beta_hat` Standard error of effect estimates
#'
#'Everything else returned has to do with simulation parameters or true effects
#'
#'+ `direct_effects`, `total_effects` matrices like G giving direct and total effects
#'
#'+ `B` true total SNP effects on traits
#'
#'+ `L_mat` true direct SNP effects on traits
#'
#'You can ignore everything else returned.
#'
#'
#'@details This function generates GWAS summary statistics from a DAG specified by the user. It can be used in two modes
#'depending on the arguments specified, 'general' mode or 'xyz' mode. Each are described below.
#'
#'In 'xyz' mode, data are generated for traits X, Y, and K variables Z_1, ..., Z_K. There is a causal effect of X on Y
#'given by gamma, which specifies the proportion of the variance of Y explained by X. The K variables Z_1, ..., Z_K can have
#'effects to or from X and Y but do not have direct effects on each other. Vectors `dir_xz` and `dir_yz` specify the direction
#'of these effects with +1 corresponding to "to" effects and -1 corresponding to "from" effects. For example, `dir_xz = c(1, -1)`
#'and `dir_yz = c(1, 1)` would indicate two variables $Z_1$ and $Z_2$ with $Z_1$ being a common cause of both $X$ and $Y$ and $Z_2$
#'being a cause of $Y$ and caused by $X$ ($Z_2$ is a mediator between $X$ and $Y$). The function will give an error if there
#'is any index with `dir_xz` equal to -1 and `dir_yz` equal 1. This would indicate a variable that is a mediator between $Y$ and
#'$X$, however, because there is an effect from $X$ to $Y$ assumed, the resulting graph would be cyclic and not allowed.
#'
#'In 'xyz' mode, you will also need to specify `tau_xz` and `tau_yz` which determine the effect sizes of each of the $Z_k$ variables on or
#'from $X$ and $Y$. These are given in signed percent variance explained. If we again used `dir_xz = c(1, -1)` and `dir_yz = c(1, 1)` with
#'`tau_xz = c(0.2, -0.3)` and `tau_yz = c(0.1, 0.25)`, this means that the confounder, $Z_1$ explains 20\% of the variance of $X$ and 10\% of the
#'variance of $Y$ and both effects are positive. $X$ explains 30\% of the variance of the mediate $Z_2$ with a negative effect direction and
#'$Z_2$ explains 25\% of the variance of $Y$. If you specify a set of variances that imply that one variable has more than 100\% of variance
#'explained you will get an error. It is also possible to specify variances that add up to less than 100\% but are incompatible with the heritabilities specified.
#'For example, if we say that $Z_1$ has a heritability of 0.8 and explains 50\% of the variance of $X$, then $X$ must have a heritability of at least 0.4.
#'
#'If using 'xyz' mode, `dir_xz`, `dir_yz`, `tau_yz`, `tau_xz` should all have the same length ($K$). `h2` and `pi` should have length equal to K+2 and have entries
#'corresponding to $X$, $Y$, and then each of the $Z_k$ variables.
#'
#'In 'general' mode, you can omit the parameters specific to 'xyz' mode and instead supply `G`, a matrix of direct effects. The diagonal
#'entries of `G` should be equal to 0 and
#'`h2` and `pi` should have length equal to the number of rows (or columns) of `G` (which should be square). All variables have variance 1
#'so `G[i,j]^2` is the proportion of variance of trait j explained by the effect of trait i.
#'Currently, this function does no checking that `G`
#'is acyclic and (I think) will run forever if `G` contains a cycle. Use caution when using this option.
#'
#'@export
sim_mv <- function(N, J,
                   tau_xz, tau_yz, dir_xz, dir_yz, gamma,
                   h2, pi, G, R_E, overlap_prop = 0,
                   R_LD = NULL, snp_info = NULL,
                   sporadic_pleiotropy = TRUE){

  if(missing(tau_xz)){
    mode <- "general"
    if(missing(G)){
      stop("If using 'general' mode, please proved G")
    }
    n <- nrow(G)
    stopifnot(ncol(G) == n)
    stopifnot(all(diag(G) == 0))
  }else{
    mode <- "xyz"
    if(missing(tau_xz) | missing(tau_yz) | missing(dir_xz)|
       missing(dir_yz) | missing(gamma)){
      stop("If using 'xyz' mode, please proved tau_xz, tau_yz, dir_xz, dir_yz, and gamma")
    }
    if(any(dir_xz == 1 & dir_yz == -1)){
      stop("No cycles allowed")
    }
    p <- length(tau_xz)
    n <- p + 2
    stopifnot(length(tau_yz) == p)
    stopifnot(length(dir_xz) == p)
    stopifnot(length(dir_yz) == p)
  }
  stopifnot(length(h2) == n | length(h2) == 1)
  if(length(h2) == 1) h2 <- rep(h2, n)
  stopifnot(length(pi) == n | length(pi) ==1)
  if(length(pi) == 1) pi <- rep(pi, n)
  stopifnot(length(N) == n | length(N) ==1)
  if(length(N) ==1) N <- rep(N, n)

  if(mode == "xyz"){
    tau_xz <- sqrt(abs(tau_xz))*sign(tau_xz)
    tau_yz <- sqrt(abs(tau_yz))*sign(tau_yz)
    gamma <- sqrt(abs(gamma))*sign(gamma)

    # Direct Effects
    G <- matrix(0, nrow = n, ncol = n)
    G[2,1] <- gamma
    G[which(dir_xz == 1) + 2, 2] <- tau_xz[dir_xz ==1]
    G[2, which(dir_xz == -1) + 2] <- tau_xz[dir_xz == -1]

    G[which(dir_yz == 1) + 2, 1] <- tau_yz[dir_yz ==1]
    G[1, which(dir_yz == -1) + 2] <- tau_yz[dir_yz == -1]
  }

  G_t <- direct_to_total(G)


  # Compute Input and Direct Heritability
  C <- colSums(h2*(G_t)^2)
  H <- t(G_t^2)
  input_h2 <- solve(diag(n) + H) %*% matrix(C, ncol = 1) %>% as.vector()

  direct_h2 <- h2-input_h2
  if(any(direct_h2 < 0)) stop("Input heritability greater total heritability for least one variable.")
  G_t <- G_t*sqrt(direct_h2)
  diag(G_t) <- sqrt(direct_h2)

  F_mat <- t(G_t)

  dat <- sim_sumstats_lf(F_mat = F_mat,
                         N = N, J = J, h_2_trait = h2,
                         omega = rep(1, n),
                         pi_L = pi,
                         overlap_prop = overlap_prop,
                         h_2_factor = rep(1, n),
                         pi_theta = 1,
                         R_E = R_E, R_LD = R_LD, snp_info  = snp_info,
                         sporadic_pleiotropy = sporadic_pleiotropy)
  direct_SNP_effects <- t(t(dat$L_mat)*diag(dat$F_mat))
  R <- list(beta_hat = dat$beta_hat,
            se_beta_hat = dat$se_beta_hat,
            direct_SNP_effects = direct_SNP_effects,
            direct_trait_effects = G,
            total_trait_effects = t(dat$F_mat)/diag(dat$F_mat),
            B = dat$Z * dat$se_beta_hat,
            R = dat$R,
            F_mat = dat$F_mat,
            dat = dat)
  if(!is.null(R_LD)) R$snp_info <- dat$snp_info

  diag(R$total_trait_effects) <- 0

  return(R)
}

path_one_step <- function(graph, paths){
  paths_new <- map(paths, function(p){
    start <- p[1]
    into <- graph[,start]
    if(all(into == 0)) return(NULL)
    new_paths <- lapply(which(into!=0), function(i){ c(i, p)})
    return(new_paths)
  })
  paths_new <- unlist(paths_new, recursive = FALSE)
  return(paths_new)
}

all_paths <- function(graph, start, end){
  p <- list(c(end))
  my_paths <- list()
  while(length(p)> 0){

    p <- path_one_step(graph, p)

    starts <- map(p, 1) %>% unlist()
    if(any(starts == start)){
      ix <- which(starts == start)
      my_paths <- append(my_paths, p[ix])
      p <- p[-ix]
    }
  }
  return(my_paths)
}

get_total_effect <- function(graph, start, end){
  paths <- all_paths(graph, start, end)
  if(length(paths) == 0) return(0)
  map(paths, function(x){
    e <- sapply(seq_along(x[-length(x)]), function(i){graph[x[i], x[i+1]]})
    return(prod(e))
  }) %>% unlist() %>%
    sum() %>%
    return()
}

direct_to_total <- function(graph){
  n = nrow(graph)
  G_total <- expand.grid(start = seq(n), end = seq(n))
  G_total$value <- map2(G_total$start, G_total$end, function(x, y){
    get_total_effect(graph, x, y)
  }) %>% unlist()
  G_total <- G_total %>% pivot_wider(names_from = end, values_from = value) %>%
    select(-start) %>% as.matrix()
  return(G_total)
}
