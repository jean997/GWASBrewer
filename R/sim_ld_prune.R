#'@title LD prune simulated data
#'@param dat Data object produced by sim_mv
#'@param pvalue A vector of p-values to use to prioritize variants.
#'@param R_LD The same list of eigen-decompositions used to generate dat
#'@param r2_thresh r^2 thershold for pruning
#'@param pval_thresh p-value threshod for pruning (see details)
#'@details Given results from sim_mv, and a vector of p-values,
#'the function will return a list of variants that have \code{p < pval_thresh}
#'and which mutually have squared correlation less than \code{r2_thresh}.
#'@return A vector of indices corresponding to the LD-pruned variant set.
#'@examples
#' data("ld_evd_list")
#' data("snpdata")
#' head(snpdata)
#' # Two traits with no causal relationship, non-overlapping GWAS
#' set.seed(1)
#' G <- matrix(0, nrow = 2, ncol = 2)
#' dat <- sim_mv(N = 10000, J = 50000, h2 = c(0.4, 0.3), pi = 1000/20000,
#'                G = G,  R_LD = ld_evd_list, snp_info = snpdata)
#'
#'# prune on p-value for first trait
#'pvals <- 2*pnorm(-abs(dat$beta_hat/dat$se_beta_hat))
#'prune_set_1 <- sim_ld_prune(dat, pvalue = pvals[,1], R_LD = ld_evd_list, pval_thresh = 1e-5)
#'length(prune_set_1)
#'@export
sim_ld_prune <- function(dat, pvalue, R_LD, r2_thresh = 0.1, pval_thresh = 1){
  if(missing(pvalue)){
    p <- 2*pnorm(-abs(dat$beta_hat/dat$se_beta_hat))
    pvalue <- apply(p, 1, min)
  }else if(length(pvalue) ==1){
    stopifnot(pvalue == round(pvalue))
    p <- 2*pnorm(-abs(dat$beta_hat/dat$se_beta_hat))
    pvalue <- p[,pvalue]
  }

  R_LD <- lapply(R_LD, function(e){e$vectors %*% diag(e$values) %*% t(e$vectors) })
  thresh <- sqrt(r2_thresh)
  combos <- expand.grid(block = unique(dat$snp_info$block),
                        rep = unique(dat$snp_info$rep))
  keep_list <- purrr::map2(combos$block, combos$rep, function(b, r){
    ix_rem <- which(dat$snp_info$block == b & dat$snp_info$rep == r)
    myld <- reshape2::melt(R_LD[[b]])
    myld$ix1 <- ix_rem[myld$Var1]
    myld$ix2 <- ix_rem[myld$Var2]
    ix_keep <- c()
    df_rem <- data.frame(ix = ix_rem, pval = pvalue[ix_rem]) %>%
      arrange(pval)
    df_rem <- filter(df_rem, pval < pval_thresh)
    while(nrow(df_rem) > 0){
      s <- df_rem$ix[1]
      ix_keep <- c(ix_keep, s)
      a <- filter(myld, (ix1 == s | ix2 == s) & abs(value) > thresh)
      v <- unique(a$ix1, a$ix2)
      df_rem <- filter(df_rem, !ix %in% v)
    }
    return(ix_keep)
  }) %>% unlist()
  return(keep_list)
}

#'@title Retrieve LD proxies
#'@param snpinfo snp_info data frame from a simulation object produced by sim_mv
#'@param index list of indexes to retrieve proxies for
#'@param R_LD list of eigen-decompositions used originally in sim_mv
#'@param r2_thresh Get proxies with r^2 >= r2_thresh with one of the index variants
#'@export
sim_ld_proxy <- function(snpinfo, index, R_LD, r2_thresh = 0.64){
  proxy_info <- lapply(index, function(i){
    b <- snpinfo$block[i]
    r <- snpinfo$rep[i]
    i <- snpinfo$ix_in_block[i]
    R <- with(R_LD[[b]], vectors %*% diag(values) %*% t(vectors))
    ii <- which(R[i,]^2 > r2_thresh)
    return(list(block = b, rep = r, index = i, proxy_index = ii, Rproxy = R[ii,ii]))
  })
  return(proxy_info)
}
