test_that("sim_mv executes", {
  # test various combinations of parameters
  # does not check that results are correct
  set.seed(98765)
  A1 <- matrix(0.7, nrow = 10, ncol = 10)
  diag(A1) <- 1
  A2 <- matrix(0.1, nrow = 6, ncol = 6)
  diag(A2) <- 1
  af <- runif(n = 16)
  # simple no LD, no AF, one study
  dat1 <- sim_mv(N = 1e5,
                 J = 100,
                 h2 = 0.05,
                 pi = 0.5,
                 G = 1)
  expect_equal(as.numeric(dat1$Sigma_G), sum(dat1$beta_joint^2))
  ## basic check on chi squared stat
  chisq <- with(dat1, sum(((beta_hat - beta_marg)/se_beta_hat)^2))
  expect_equal(chisq, 83.181452) # this will fail if anything updates but then i will know to re-run checks
  # expect_lt(chisq, qchisq(p = 0.005, df = 100, lower.tail = F))
  # expect_gt(chisq, qchisq(p = 0.005, df = 100, lower.tail = T))
  expect_equal(dat1$beta_joint, dat1$beta_marg)
  expect_equal(dat1$beta_joint, dat1$direct_SNP_effects_joint)
  expect_equal(dat1$beta_joint, dat1$direct_SNP_effects_marg)
  expect_equal(as.numeric(dat1$Sigma_G + dat1$Sigma_E), 1)
  ## check pi_exact and h2_exact
  dat2 <- sim_mv(N = 1e5,
                 J = 100,
                 h2 = 0.05,
                 pi = 0.5,
                 G = 1,
                 pi_exact = TRUE,
                 h2_exact = TRUE)
  expect_equal(sum(dat2$beta_joint != 0), 50)
  expect_equal(as.numeric(dat2$Sigma_G), 0.05)

  ## check sample size specifications and R_E/R_obs
  N <- matrix(c(4e5, 1e5, 1e5, 2e5), nrow = 2)
  Nc <- diag(2)
  Nc[1,2] <- Nc[2,1] <- 1/(sqrt(2)*sqrt(4))
  R_E <- matrix(0.8, nrow = 2, ncol = 2)
  diag(R_E) <- 1
  dat3 <- sim_mv(N = N,
                 J = 100,
                 h2 = 0.05,
                 pi = 0.5,
                 G = 2,
                 R_E = R_E)
  expect_equal(stats::cov2cor(dat3$Sigma_E), R_E)
  expect_equal(dat3$R, Nc*dat3$trait_corr)

  set.seed(5)
  dat4 <- sim_mv(N = N,
                 J = 100,
                 h2 = 0.05,
                 pi = 0.5,
                 G = 2,
                 R_obs = R_E)
  expect_equal(dat4$trait_corr, R_E)
  expect_equal(dat4$R, Nc*dat4$trait_corr)
  # data frame sample size specification
  N_df <- data.frame("trait_1" = c(1, 1, 0), "trait_2" = c(0, 1, 1), "N" = c(3e5, 1e5, 1e5))
  set.seed(5)
  dat5 <- sim_mv(N = N_df,
                 J = 100,
                 h2 = 0.05,
                 pi = 0.5,
                 G = 2,
                 R_obs = R_E)
  expect_identical(dat4, dat5)
  # zero sample size for missing sumstats
  N[2,] <- 0
  N[,2] <- 0
  dat6 <- sim_mv(N = N,
                 J = 100,
                 h2 = 0.05,
                 pi = 0.5,
                 G = 2,
                 R_obs = R_E)
  expect_equal(dat6$beta_hat[,2], rep(NA_real_, 100))
  expect_equal(dat6$se_beta_hat[,2], rep(NA_real_, 100))

  ## simple LD pattern
  set.seed(5566)
  dat7 <- sim_mv(N = N,
                 J = 20,
                 h2 = 0.05,
                 pi = 1,
                 G = 2,
                 R_obs = R_E,
                 R_LD = list(A1, A2),
                 af = af)
  expect_equal(dat7$snp_info$AF, rep(af, 2)[1:20])
  expect_equal(dat7$Sigma_G + dat7$Sigma_E, dat7$trait_corr)
  L <- as.matrix(Matrix::bdiag(A1, A2, A1[1:4, 1:4]))
  sx <- with(dat7$snp_info, sqrt(2*AF*(1-AF)))
  S <- diag(1/sx)
  expect_equal(dat7$beta_marg, S %*% L %*% solve(S) %*% dat7$beta_joint)
  expect_equal(dat7$direct_SNP_effects_marg, S %*% L %*% solve(S) %*% dat7$direct_SNP_effects_joint)
  expect_equal(dat7$beta_joint, dat7$direct_SNP_effects_joint)
  expect_equal(dat7$beta_marg, dat7$direct_SNP_effects_marg)
  expect_equal(dat7$beta_hat[1,1], 0.1, tolerance = 1e-1) # this just tells me if gen_bhat_from_b has changed
})
