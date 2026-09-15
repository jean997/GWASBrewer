#' Example allele frequencies
#'
#' A vector of allele frequencies for variants on chromosome 19,
#' provided for use in examples and vignettes.
#'
#' @format A numeric vector of length 19490.
"AF"

#' Example LD block matrices
#'
#' A list of LD block correlation matrices used to demonstrate
#' simulation with LD.
#'
#' @format A list of 39 sparse correlation matrices (class \code{dsCMatrix}).
"ld_mat_list"

#' Example SNP data
#'
#' Variant information accompanying \code{\link{AF}} and \code{\link{ld_mat_list}},
#' provided for use in examples and vignettes.
#'
#' @format A data frame with 19490 rows and 14 variables including allele frequency,
#' rs ID, chromosome, position, LD block id, and LD pruning indicators.
"snpdata"
