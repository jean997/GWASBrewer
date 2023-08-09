#'@title LD prune simulated data
#'@param dat Data object produced by sim_mv
#'@param pvalue Either a vector used to prioritize variants or an integer.
#'If pvalue is an integer, i, variants will be priorized by the p-value for trait i.
#'If pvalue is missing, variants will be prioritized randomly.
#'@param R_LD LD pattern used to generate dat
#'@param r2_thresh r^2 threshold for pruning
#'@param pval_thresh p-value threshold for pruning (see details)
#'@details Given results from sim_mv, and a vector of p-values,
#'the function will return a list of variants that have \code{p < pval_thresh}
#'and which mutually have squared correlation less than \code{r2_thresh}.
#'@return A vector of indices corresponding to the LD-pruned variant set.
#'@examples
#' data("ld_mat_list")
#' data("snpdata")
#'
#' # Two traits with no causal relationship, non-overlapping GWAS
#' set.seed(1)
#' G <- matrix(0, nrow = 2, ncol = 2)
#' dat <- sim_mv(N = 10000, J = 50000, h2 = c(0.4, 0.3), pi = 1000/20000,
#'                G = G,  R_LD = ld_mat_list, snp_info = snpdata)
#'
#'# prune on p-value for first trait
#'pvals <- 2*pnorm(-abs(dat$beta_hat/dat$se_beta_hat))
#'prune_set_1 <- sim_ld_prune(dat, pvalue = pvals[,1], R_LD = ld_mat_list, pval_thresh = 1e-5)
#'# Above is equivalent to
#'prune_set_1 <- sim_ld_prune(dat, pvalue = 1, R_LD = ld_mat_list, pval_thresh = 1e-5)
#'@export
sim_ld_prune <- function(dat, pvalue, R_LD, r2_thresh = 0.1, pval_thresh = 1){
  stopifnot("sim_mv" %in% class(dat))
  if(missing(pvalue)){
    if(!pval_thresh == 1){stop("If providing pval_thresh, please also provide pvalue.\n")}
    message("pvalue omitted so variants will be prioritized randomly")
    pvalue <- runif(n = nrow(dat$beta_hat))
  }else if(length(pvalue) ==1){
    stopifnot(pvalue == round(pvalue))
    message(paste0("Prioritizing variants based on p-value for trait ", pvalue, "\n"))
    if(!is.null(dat$s_estimate)){
      p <- 2*pnorm(-abs(dat$beta_hat/dat$s_estimate))
    }else{
      p <- 2*pnorm(-abs(dat$beta_hat/dat$se_beta_hat))
    }
    pvalue <- p[,pvalue]
  }


  R_LD <- check_R_LD(R_LD, return = "matrix")
  ii <- lapply(R_LD, function(R){
    which(R^2 > r2_thresh, arr.ind = TRUE) %>%
      data.frame()
  })

  thresh <- sqrt(r2_thresh)
  combos <- expand.grid(block = unique(dat$snp_info$block),
                        rep = unique(dat$snp_info$rep))
  dat$snp_info$ix_in_dat <- seq(nrow(dat$snp_info))

  keep_list <- purrr::map2(combos$block, combos$rep, function(b, r){
    d <- filter(dat$snp_info, block == b & rep == r) %>% select(ix_in_block, ix_in_dat) %>% arrange(ix_in_block)
    d$pval <- pvalue[d$ix_in_dat]
    df_rem <- d %>%
      arrange(pval) %>%
      filter(pval < pval_thresh)
    # myld <- melt_to_df(R_LD[[b]]^2)
    # myld$ix_in_dat1 <- d$ix_in_dat[myld$row]
    # myld$ix_in_dat2 <- d$ix_in_dat[myld$col]
    # myld <- filter(myld, ix_in_dat1 %in% df_rem$ix_in_dat & ix_in_dat2 %in% df_rem$ix_in_dat)
    # ix_keep <- c()
    # while(nrow(df_rem) > 0){
    #   s <- df_rem$ix_in_dat[1]
    #   ix_keep <- c(ix_keep, s)
    #   a <- filter(myld, (ix_in_dat1 == s | ix_in_dat2 == s) & abs(value) > thresh)
    #   v <- unique(a$ix_in_dat1, a$ix_in_dat2)
    #   df_rem <- filter(df_rem, !ix_in_dat %in% v)
    # }

    ix_keep <- c()
    while(nrow(df_rem) > 0){
      s <- df_rem$ix_in_dat[1]
      sb <- df_rem$ix_in_block[1]
      ix_keep <- c(ix_keep, s)
      a <- ii[[b]]$col[ii[[b]]$row == sb]
      df_rem <- filter(df_rem, !ix_in_block %in% a)
    }

    return(ix_keep)
  }) %>% unlist()
  return(keep_list)
}

# A faster version of melt
melt_to_df <- function(x){
  data.frame(row=seq(nrow(x))[row(x)],
             col=seq(ncol(x))[col(x)],
             value=as.numeric(x))
}

#'@title Retrieve LD proxies
#'@param dat Simulation object produced by `sim_mv`
#'@param index list of indexes to retrieve proxies for
#'@param R_LD LD pattern used to generate dat
#'@param r2_thresh Get proxies with r^2 >= r2_thresh with one of the index variants
#'@param return_mat If TRUE, return the correlation matrix between the index variant and the proxies. In this matrix, the
#'the index variant always corresponds the first row/column and the proxies are in the order returned.
#'@export
sim_ld_proxy <- function(dat, index, R_LD, r2_thresh = 0.64, return_mat = FALSE){
  R_LD <- check_R_LD(R_LD, "matrix")
  dat$snp_info$ix_in_dat <- seq(nrow(dat$snp_info))
  ii <- list()
  blocks <- unique(dat$snp_info[index,]$block)
  ii[blocks] <- lapply(blocks, function(b){
    which(R_LD[[b]]^2 > r2_thresh, arr.ind = TRUE) %>%
      data.frame() %>%
      filter(row  != col)
  })
  proxy_info <- lapply(index, function(i){
    b <- dat$snp_info$block[i]
    r <- dat$snp_info$rep[i]
    id <- dat$snp_info$ix_in_block[i]
    proxy_index_block <- ii[[b]]$col[ii[[b]]$row == id] %>% sort()
    proxy_index_dat <- filter(dat$snp_info, rep == r & block == b & ix_in_block %in% proxy_index_block) #%>%
                       #arrange(ix_in_block) %>%
                       #with(., ix_in_dat)
    ret <- list( index = i,
                # block = b, rep = r,
                #proxy_index_block = proxy_index_block,
                proxy_index = proxy_index_dat$ix_in_dat)
    if(return_mat){
      if (length(ret$proxy_index) == 0) {
          print("no proxy SNPs found!")
          ret$Rproxy <- matrix(1, nrow = 1)
        }
      else {
          ret$Rproxy <- R_LD[[b]][c(id,proxy_index_dat$ix_in_block), c(id,proxy_index_dat$ix_in_block)]
          rownames(ret$Rproxy) <- colnames(ret$Rproxy) <-  c(i,proxy_index_dat$ix_in_dat)
        }
    }
    return(ret)
  })
  return(proxy_info)
}

#'@title Extract LD matrix from simulated data
#'@description Extract LD matrix for specific variants in simulated data set
#'@param dat Simulation object produced by `sim_mv`
#'@param index vector of indices for snps to extract LD
#'@param R_LD List of eigen-decompositions used in original simulation
#'@return An LD matrix. SNP order matches original index order
#'@examples
#' data("ld_mat_list")
#' data("snpdata")
#' # Two traits with no causal relationship, non-overlapping GWAS
#' set.seed(1)
#' G <- matrix(0, nrow = 2, ncol = 2)
#' dat <- sim_mv(N = 10000, J = 50000, h2 = c(0.4, 0.3), pi = 1000/20000,
#'                G = G,  R_LD = ld_mat_list, snp_info = snpdata)
#'
#'# extract ld matrix for all variants with p-value for trait 1 less than 1e-5
#'pvals <- 2*pnorm(-abs(dat$beta_hat/dat$se_beta_hat))
#'index <- which(pvals[,1] < 1e-5)
#'ld_mat <- sim_extract_ld(dat, index, ld_mat_list)
#'dim(ld_mat)
#'length(index)
#'@export
sim_extract_ld <- function(dat, index, R_LD){
  if(!length(index) == length(unique(index))){
    stop("index should contain only unique indices (no duplicates)\n")
  }
  R_LD <- check_R_LD(R_LD, "matrix")
  snpinfo <- dat$snp_info[index,]
  block_sum <- snpinfo %>% group_by(block, rep) %>% summarize(n = n())
  ld_sub_blocks <- map(unique(block_sum$block), function(b){
    r <- block_sum$rep[block_sum$block == b]
    sub_sub_blocks <- map(r, function(r){
      i <- filter(snpinfo, block == b& rep == r) %>% with(., ix_in_block)
      s <- filter(snpinfo, block == b& rep == r) %>% with(., SNP)
      return(list("mat" = R_LD[[b]][i, i, drop = FALSE], "snp" = s))
    })
    return(sub_sub_blocks)
  }) %>% unlist(recursive = FALSE)
  mats <- map(ld_sub_blocks, "mat")
  snps <- map(ld_sub_blocks, "snp") %>% unlist()
  mat <- Matrix::bdiag(mats) %>% as.matrix()
  o <- match(snpinfo$SNP, snps)
  snps <- snps[o]
  mat <- mat[o, o]
  rownames(mat) <- colnames(mat) <- index
  return(mat)
}
