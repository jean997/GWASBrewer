#'@title Simulate multivariate GWAS data
#'@param N GWAS sample size. N can be a scalar, vector, or matrix. If N is a scalar, all GWAS have the same sample size
#'and there is no overlap between studies. If N is a vector, each element of N specifies the sample size of the corresponding
#'GWAS and there is no overlap between studies. If N is a matrix, N_ii specifies the sample size of study i
#' and N_ij specifies the number of samples present in both study i and study j. The elements of N must be positive
#' but non-integer values will not generate an error.
#'@param J Number of variants to simulate
#'@param h2 A scalar or vector giving the heritability of each trait.
#'@param pi A scalar or vector giving the expected proportion of direct effect SNPs for each trait.
#'@param G Matrix of direct effects. Rows correspond to the 'from' trait
#'and columns correspond to the 'to' trait, so \code{G[1,2]} is the direct effect of trait 1 on trait 2. G should have 0
#'on the diagonal.
#'@param R_E Environmental correlation between traits. R_E is ignored if there is no sample overlap.
#'@param R_LD List of LD blocks. R_LD should have class \code{list}.
#'Each element of R_LD can be either a) a matrix, b) a sparse matrix (class \code{dsCMatrix}) or c) an eigen decomposition (class \code{eigen}).
#'All elements be correlation matrices, meaning that they have 1 on the diagonal and are positive definite. See Details.
#'@param af Optional vector of allele frequencies. If R_LD is not supplied, af can be a scalar, vector or function.
#'If af is a function it should take a single argument (n) and return a vector of n allele frequencies (See Examples).
#'If R_LD is supplied, af must be a vector with length equal to the size of the supplied LD pattern (See Examples).
#'@param sporadic_pleiotropy Allow sporadic pleiotropy between traits. Defaults to TRUE.
#'@param pi_exact If TRUE, the number of direct effect SNPs for each trait will be exactly equal to round(pi*J).
#'@param h2_exact If TRUE, the heritability of each trait will be exactly h2.
#'@param est_s If TRUE, return estimates of se(beta_hat).
#'@param return_dat Useful development option.
#'
#'@return A list with the following elements:
#'
#'Simulated effect estimates and standard errors are contained in matrices
#'
#' + \code{beta_hat}: Effect estimates for each trait
#'
#' + \code{se_beta_hat} Standard error of effect estimates, equal to sqrt(1/N_m *Var(G_j)).
#'
#' + \code{s_estimate} Estimate of \code{se_beta_hat}. Present only if \code{est_s = TRUE}.
#'
#' Four matrices contain direct and total marginal and joint SNP-trait associations:
#' + \code{direct_SNP_effects_marg} and \code{direct_SNP_effects_joint} give direct effects of SNPs on traits. These are the
#'same if there is no LD. If there is LD, \code{direct_SNP_effects_marg} is the direct component of the expected marginal association.
#'
#' + \code{beta_marg} and \code{beta_joint} give the total SNP effects (direct and indirect). These are the
#'same if there is no LD. If there is LD, \code{beta_marg} is the total expected marginal association. i.e. \code{beta_marg} is the
#'expected value of \code{beta_hat}.
#'
#' Four matrices describe the covariance of traits and the row correlation of effect estimates
#'
#' + \code{Sigma_G} Genetic variance-covariance matrix. Diagonal elements are equal to heritability.
#'
#' + \code{Sigma_E} Environmental variance-covariance matrix.
#'
#' + \code{trait_corr} Population correlation of traits. \code{trait_corr = Sigma_G + Sigma_E}.
#'
#' + \code{R} Row correlation of \code{beta_hat - beta_marg}, equal to \code{trait_corr} scaled by the overlap proportion matrix.
#'
#' Finally,
#'
#' + \code{snp_info} Is a data frame with variant information. If \code{R_LD} was omitted, \code{snp_info} contains only the allele frequency of
#' each variant. If \code{R_LD} was included, \code{snp_info} also contains block and replicate information.
#'
#'@details This function generates GWAS summary statistics from a linear SEM specified by the matrix G.
#'The previous "xyz" mode is now deprecated. If you are used to using this function in xyz mode,
#'you can generate the corresponding G using the xyz_to_G function (see examples).
#'
#'G should be square (#traits by #traits) matrix with 0s on the diagonal.
#'All traits have variance 1, so \code{G[i,j]^2} is the proportion of variance of trait j explained by the effect of trait i.
#'You will get an error if you specify a cyclic DAG, a DAG that is impossible, or a combination of \code{G} and \code{h2}
#'that is impossible.
#'
#'To generate data with LD, supply \code{R_LD} which can be a list of matrices, sparse matrices, or eigen-decompositions.
#'These matrices are interpreted as blocks in a block-diagonal LD matrix. The \code{af} argument must be provided
#'and should be a vector with length equal to the total size of the LD pattern (the sum of the sizes of each block).
#'\code{R_LD} does not need to have the same size as \code{J}. The LD pattern will be repeated or subset as necessary
#'to generate the desired number of variants.
#'
#'@examples
#' # Two traits with no causal relationship and some environmental correlation
#' # specify completely overlapping GWAS
#' N <- matrix(1000, nrow = 2, ncol =2)
#' G <- matrix(0, nrow = 2, ncol = 2)
#' R_E <- matrix(c(1, 0.8, 0.8, 1), nrow = 2, ncol = 2)
#' dat <- sim_mv(N = N, J = 20000, h2 = c(0.4, 0.3), pi = 1000/20000,
#'                G = G, R_E = R_E)
#'dat$R
#'cor(dat$beta_hat - dat$beta_marg)
#'
#'# The af argument can be a scalar, vector, or function.
#'dat <- sim_mv(N = N, J = 20000, h2 = c(0.4, 0.3), pi = 1000/20000,
#'                G = G, R_E = R_E, af = function(n){rbeta(n = n, 1, 5)})
#'
#'# A very simple example with LD
#'# Use a pattern of two small blocks of LD
#'A1 <- matrix(0.7, nrow = 10, ncol = 10)
#'diag(A1) <- 1
#'A2 <- matrix(0.1, nrow = 6, ncol = 6)
#'diag(A2) <- 1
#'# If using LD, af should have the same size as the LD pattern
#'af <- runif(n = 16)
#'dat <- sim_mv(N = N, J = 20000, h2 = c(0.4, 0.3), pi = 1000/20000,
#'                G = G, R_E = R_E, R_LD = list(A1, A2), af = af)
#'
#' # Use xyz_to_G to generate G from xyz specification
#' myG <- xyz_to_G(tau_xz = c(0.2, -0.3), tau_yz = c(0.1, 0.25),
#'         dir_xz = c(1, -1), dir_yz = c(1,1), gamma = 0)
#' # If N is a scalar or a vector, there is no sample overlap
#' dat <- sim_mv(N = 10000, J = 20000, h2 = rep(0.4, 4),
#'               pi = c(500, 500, 1000, 1000)/20000,
#'               G = myG)
#'plot(dat$beta_marg[,3], dat$beta_marg[,1])
#'abline(0, dat$total_trait_effects[3,1])
#'@export
sim_mv <- function(N, J,
                   h2, pi, G=0, R_E = NULL,
                   R_LD = NULL, af = NULL,
                   snp_effect_function = "normal",
                   sporadic_pleiotropy = TRUE,
                   pi_exact = FALSE,
                   h2_exact = FALSE,
                   est_s = FALSE,
                   return_dat  = FALSE){


  if(!"matrix" %in% class(G)) n <- 1
    else n <- nrow(G)
  G <- check_matrix(G, "G", n, n)
  h2 <- check_scalar_or_numeric(h2, "h2", n)
  h2 <- check_01(h2)
  pi <- check_pi(pi, J, n)
  nn <- check_N(N, n)
  G <- check_G(G, h2, n)

  G_t <- G$G_tot*sqrt(G$dir_h2)
  diag(G_t) <- sqrt(G$dir_h2)

  F_mat <- t(G_t)

  dat <- sim_lf(F_mat = F_mat,
                N = N, J = J,
                h2_trait = h2,
                omega = rep(1, n),
                h2_factor = rep(1, n),
                pi_L = pi, pi_theta = 1,
                R_LD = R_LD, af = af,
                R_E = R_E,
                snp_effect_function = snp_effect_function,
                sporadic_pleiotropy = sporadic_pleiotropy,
                est_s = est_s,
                h2_exact = h2_exact,
                pi_exact = pi_exact)

  direct_SNP_effects <- t(t(dat$L_mat)*diag(dat$F_mat))

  direct_SNP_effects_joint <- t(t(dat$L_mat_joint)*diag(dat$F_mat))
  R <- list(beta_hat = dat$beta_hat,
            se_beta_hat = dat$se_beta_hat,
            direct_SNP_effects_marg = direct_SNP_effects,
            direct_SNP_effects_joint = direct_SNP_effects_joint,
            direct_trait_effects = G$G_dir,
            total_trait_effects = G$G_tot,
            beta_joint = dat$beta_joint,
            beta_marg = dat$beta_marg,
            trait_corr = dat$trait_corr,
            R = dat$R,
            Sigma_G = dat$Sigma_G,
            Sigma_E = dat$Sigma_E,
            snp_info = dat$snp_info)
  if(est_s){
    R$s_estimate <- dat$s_estimate
  }
  if(return_dat){
    R$L_mat <- dat$L_mat
    R$L_mat_joint_std <- dat$L_mat_joint_std
    R$F_mat <- dat$F_mat
  }

  return(R)
}


