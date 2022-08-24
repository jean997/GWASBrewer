#'@export
sim_ld_prune <- function(dat, pvalue, R_LD, r2_thresh = 0.1, pval_thresh = 1){
  if(missing(pvalue)){
    p <- 2*pnorm(-abs(dat$beta_hat/dat$se_beta_hat))
    pvalue <- apply(p, 1, min)
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
