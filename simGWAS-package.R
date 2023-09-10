#' simGWAS: Simulate GWAS summary statistics from specified DAG or factor structure.
#' 
#' @importFrom stats cov2cor pnorm qnorm rbinom rnorm runif cov2cor
#' @importFrom utils str
#' @docType package
#' @name simGWAS
NULL

#' SNP data for examples
#' 
#' The package contains a built-in data set containing the LD pattern from Chromosome 19
#' in HapMap3 broken into 39 blocks. This LD pattern was estimated from the HapMap3 European
#' subset using LDShrink. This data set can also be downloaded 
#' [here](https://zenodo.org/record/6761943#.Yrno2njMIUE). The LD pattern must be accompanied
#' by a vector of allele frequencies with length equal to the total size of the LDpattern 
#' (i.e. the sum of the size of each block in the list).
#' 
#' @format ## `snpdata`
#' A tibble with 19490 rows and 14 columns:
#' \describe{
#'   \item{AF}{AF}
#'   \item{SNP}{SNP}
#'   \item{allele}{allele}
#'   \item{chr}{chr}
#'   \item{ld_snp_id}{ld_snp_id}
#'   \item{map}{map}
#'   \item{pos}{pos}
#'   \item{region_id}{region_id}
#'   \item{snp_id}{snp_id}
#'   \item{in_hapmap}{in_hapmap}
#'   \item{ldscore_1kg}{ldscore_1kg}
#'   \item{ldscore_hm3}{ldscore_hm3}
#'   \item{keep_ld_prune_0}{keep_ld_prune_0}
#'   \item{keep_ld_prune_0}{keep_ld_prune_0}
#' }
#' @source <https://zenodo.org/record/6761943#.Yrno2njMIUE>
"snpdata"

#' Allele frequencies for SNP data
#' 
#' Allele frequency column from `snpdata` object
"AF"

#' Block LD matrices for SNP data
"ld_mat_list"
