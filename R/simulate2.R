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
#'and columns correspond to the 'to' trait, so G[1,2] is the direct effect of trait 1 on trait 2. G should have 0
#'on the diagonal.
#'@param R_E Environmental correlation between traits. R_E is ignored if there is no sample overlap.
#'@param R_LD List of eigen decompositions of LD correlation matrices, may be missing.
#'@param snp_info If R_LD is provided, provide a data frame with columns "SNP" and "AF"
#'@param af Optional vector of allele frequencies.
#'This option is only used when data are generated without LD.
#'Otherwise the allele frequency in snp_info is used instead.
#'@param sporadic_pleiotropy Allow sporadic pleiotropy between traits. Defaults to TRUE.
#'@param return_dat Useful for debugging.
#'
#'@return A list with the following elements:
#'
#'Simulated effect estimates and standard errors are contained in matrices
#'
#' + \code{beta_hat}: Effect estimates for each trait
#'
#' + \code{se_beta_hat} Standard error of effect estimates
#'
#'Everything else returned has to do with simulation parameters or true effects
#'
#'+ \code{direct_SNP_effects_marg} and \code{direct_SNP_effects_joint} give direct effects of SNPs on traits. These are the
#'same if there is no LD. If there is LD, \code{direct_SNP_effects_marg} is the direct component of the expected marginal association.
#'
#'+ \code{beta_marg} and \code{beta_joint} give the total SNP effects (direct and indirect). These are the
#'same if there is no LD. If there is LD, \code{beta_marg} is the total expected marginal association. i.e. \code{beta_marg} is the
#'expected value of \code{beta_hat}.
#'
#'+ \code{direct_trait_effects} and \code{total_trait_effects} matrices giving direct and total effects of traits on each other.
#'\code{direct_trait_effects} should be equal to the supplied \code{G}.
#'
#'+ \code{R} gives the row correlation of \code{beta_hat - beta_marg}.
#'
#'Most users can ignore everything else returned.
#'
#'
#'@details This function generates GWAS summary statistics from a linear SEM specified by the matrix G.
#'The previous "xyz" mode is now deprecated. If you are used to using this function in xyz mode,
#'you can generate the corresponding G using the xyz_to_G function (see examples).
#'
#'G should be square nxn matrix with 0s on the diagonal. \code{h2} and \code{pi} should be scalars or have length n.
#'All variables have variance 1, so \code{G[i,j]^2} is the proportion of variance of trait j explained by the effect of trait i.
#'
#'You should get an error if you specify a cyclic DAG, a DAG that is impossible, or a combination of \code{G} and \code{h2}
#'that is impossible.
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
#'# A very simple example with LD
#'# Use a pattern of two small blocks of LD
#'A1 <- matrix(0.7, nrow = 10, ncol = 10)
#'diag(A1) <- 1
#'A2 <- matrix(0.1, nrow = 6, ncol = 6)
#'diag(A2) <- 1
#'snp_info <- data.frame(SNP = 1:16, AF = runif(n = 16))
#'dat <- sim_mv(N = N, J = 20000, h2 = c(0.4, 0.3), pi = 1000/20000,
#'                G = G, R_E = R_E, R_LD = list(A1, A2), snp_info = snp_info)
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
                   h2, pi, G, R_E = NULL,
                   R_LD = NULL, snp_info = NULL, af = NULL,
                   sporadic_pleiotropy = TRUE,
                   estimate_s = FALSE,
                   return_dat  = FALSE){


  n <- nrow(G)
  G <- check_matrix(G, "G", n, n)
  h2 <- check_scalar_or_numeric(h2, "h2", n)
  h2 <- check_01(h2)
  pi <- check_scalar_or_numeric(pi, "pi", n)
  nn <- check_N(N, n)
  G <- check_G(G, h2, n)

  G_t <- G$G_tot*sqrt(G$dir_h2)
  diag(G_t) <- sqrt(G$dir_h2)

  F_mat <- t(G_t)

  dat <- sim_sumstats_lf(F_mat = F_mat,
                         N = N, J = J,
                         h2_trait = h2,
                         omega = rep(1, n),
                         pi_L = pi,
                         af = af,
                         h2_factor = rep(1, n),
                         pi_theta = 1,
                         R_E = R_E,
                         R_LD = R_LD,
                         snp_info  = snp_info,
                         sporadic_pleiotropy = sporadic_pleiotropy,
                         estimate_s = estimate_s)
  direct_SNP_effects <- t(t(dat$L_mat)*diag(dat$F_mat))

  direct_SNP_effects_joint <- t(t(dat$L_mat_joint)*diag(dat$F_mat))
  R <- list(beta_hat = dat$beta_hat,
            se_beta_hat = dat$se_beta_hat,
            direct_SNP_effects_marg = direct_SNP_effects,
            direct_SNP_effects_joint = direct_SNP_effects_joint,
            direct_trait_effects = G$G_dir,
            total_trait_effects = t(dat$F_mat)/diag(dat$F_mat),
            beta_joint = dat$beta_joint,
            beta_marg = dat$beta_marg,
            trait_corr = dat$trait_corr,
            R = dat$R,
            R_E = dat$R_E,
            true_h2 = dat$true_h2,
            af = dat$af)
  if(estimate_s){
    R$s_estimate <- dat$s_estimate
  }
  if(return_dat){
    R$L_mat <- dat$L_mat
    R$L_mat_joint_std <- dat$L_mat_joint_std
    R$F_mat <- dat$F_mat
  }
  if(!is.null(R_LD)) R$snp_info <- dat$snp_info
  diag(R$total_trait_effects) <- 0
  return(R)
}


