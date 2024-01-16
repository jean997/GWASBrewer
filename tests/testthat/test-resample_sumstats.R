test_that("resample_sumstats executes", {
  set.seed(98765)
  A1 <- matrix(0.7, nrow = 10, ncol = 10)
  diag(A1) <- 1
  A2 <- matrix(0.1, nrow = 6, ncol = 6)
  diag(A2) <- 1
  af1 <- runif(n = 16, min = 0.3, max = 0.7)
  af2 <- runif(n = 16, min = 0.05, max = 0.25)
  af <- stats::rbeta(n = 100, 1, 5)
  # simple no LD, no AF, one study
  dat1 <- sim_mv(N = 1e5,
                 J = 100,
                 h2 = 0.05,
                 pi = 0.5,
                 G = 2)
  dat1_1 <- resample_sumstats(dat = dat1,
                              N = 3e5,
                              geno_scale = "sd")
  expect_identical(dat1$beta_joint, dat1_1$beta_joint)
  expect_identical(dat1$beta_marg, dat1_1$beta_marg)
  #chisq <- with(dat1_1, sum(((beta_hat - beta_marg)/se_beta_hat)^2))

  dat1_rescale <- rescale_sumstats(dat = dat1, output_geno_scale = "allele", af = af)
  sx <- sqrt(2*af*(1-af))
  expect_equal(dat1$beta_joint, dat1_rescale$beta_joint*sx)
  set.seed(1)
  dat1_2 <- resample_sumstats(dat = dat1, N  = 3e5, af = af, geno_scale = "allele")
  set.seed(1)
  dat1_3 <- resample_sumstats(dat = dat1_rescale, N = 3e5, af = af, geno_scale = "allele")
  expect_identical(dat1_2, dat1_3, tolerance = 1e-15)

  dat1_rescale <- rescale_sumstats(dat = dat1, output_geno_scale = "allele", af = af, output_pheno_sd = 1.4)
  expect_equal(diag(dat1_rescale$Sigma_G + dat1_rescale$Sigma_E), rep(1.4^2, 2))
  set.seed(1)
  dat1_4 <- resample_sumstats(dat = dat1_rescale, N = 3e5, af = af, geno_scale = "allele")
  dat1_4 <- rescale_sumstats(dat = dat1_4, output_geno_scale = "allele", output_pheno_sd = 1)
  expect_identical(dat1_2, dat1_4, tolerance = 1e-15)

  set.seed(1)
  dat1_5 <- resample_sumstats(dat = dat1, N = 3e5, af = af,
                              geno_scale = "allele",
                              new_env_var = c(0.7, 1.3))
  expect_equal(dat1_5$h2, diag(dat1$Sigma_G)/dat1_5$pheno_sd^2)
  expect_equal(diag(dat1_5$Sigma_G) + diag(dat1_5$Sigma_E), dat1_5$pheno_sd^2)
  expect_equal(dat1_5$pheno_sd^2, c(0.7, 1.3) + diag(dat1$Sigma_G))
  expect_equal(dat1_2$beta_joint, dat1_5$beta_joint)
  expect_equal(dat1_2$se_beta_hat, GWASBrewer:::row_times(dat1_5$se_beta_hat, 1/dat1_5$pheno_sd))

  ## allele frequncies and LD different in two populations
  set.seed(1000)
  G <- matrix(0, nrow = 2, ncol = 2)
  G[1,2] <- 0.8
  dat2 <- sim_mv(N = 1e5,
                 J = 17,
                 h2 = 0.2,
                 pi = 1,
                 G = G,
                 R_LD = list(A1, A2),
                 af = af1)
  expect_equal(dat2$beta_joint[,1]*0.8 + dat2$direct_SNP_effects_joint[,2], dat2$beta_joint[,2])
  L <- as.matrix(Matrix::bdiag(A1, A2, matrix(1)))
  sx <- sqrt(2*af1*(1-af1))[c(1:16, 1)]
  S <- diag(1/sx)
  expect_equal(dat2$beta_marg, S %*% L %*% solve(S) %*% dat2$beta_joint)
  expect_equal(dat2$direct_SNP_effects_marg, S %*% L %*% solve(S) %*% dat2$direct_SNP_effects_joint)


  R_obsn <- matrix(-0.8, nrow = 2, ncol = 2)
  diag(R_obsn) <- 1
  expect_error(dat2_1 <- resample_sumstats(dat = dat2,
                                           R_LD = list(A1, A2),
                                           af = af1,
                                           new_R_obs = R_obsn,
                                           N = 3e5))

  dat2_1 <- resample_sumstats(dat = dat2, R_LD = list(A1, A2), af = af1, N = 3e5)
  expect_equal(dat2$beta_joint, dat2_1$beta_joint)
  expect_equal(dat2$beta_marg, dat2_1$beta_marg)
  expect_equal(dat2$direct_SNP_effects_marg, dat2_1$direct_SNP_effects_marg)

  dat2_2 <- resample_sumstats(dat = dat2, R_LD = list(A2, A1), af = af2, N = 3e5)
  expect_equal(dat2$beta_joint, dat2_2$beta_joint)
  expect_equal(dat2$Sigma_E, dat2_2$Sigma_E)
  expect_equal(dat2_2$h2, diag(dat2_2$Sigma_G)/dat2_2$pheno_sd^2)
  expect_equal(diag(dat2_2$Sigma_G) + diag(dat2_2$Sigma_E), dat2_2$pheno_sd^2)
  expect_lt(dat2_2$pheno_sd[1], dat2$pheno_sd[1])
  L <- as.matrix(Matrix::bdiag(A2, A1, matrix(1)))
  sx <- sqrt(2*af2*(1-af2))[c(1:16, 1)]
  S <- diag(1/sx)
  expect_equal(dat2_2$beta_marg, S %*% L %*% solve(S) %*% dat2_2$beta_joint)
  expect_equal(dat2_2$direct_SNP_effects_marg, S %*% L %*% solve(S) %*% dat2_2$direct_SNP_effects_joint)
})
