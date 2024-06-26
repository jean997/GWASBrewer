% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_h2.R
\name{compute_h2}
\alias{compute_h2}
\title{Compute heritability from standardized or non-standardized effects}
\usage{
compute_h2(
  b_joint,
  geno_scale = c("allele", "sd"),
  pheno_sd = 1,
  R_LD = NULL,
  af = NULL,
  full_mat = FALSE
)
}
\arguments{
\item{R_LD}{LD pattern (optional). See \code{?sim_mv} for more details.}

\item{af}{Allele frequencies (optional, allowed only if \code{R_LD} is missing). See \code{?sim_mv} for more details.}

\item{full_mat}{If TRUE, return the full genetic variance-covariance matrix}

\item{b_joint_std, b_joint}{matrix of standardized or non-standardized effects. Provide only one of these options.}
}
\description{
Compute heritability from standardized or non-standardized effects
}
