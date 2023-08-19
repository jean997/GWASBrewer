test_that("resample_inddata executes", {
  set.seed(17576)
  A <- matrix(c( 1.00,  0.60,  0.40, -0.1, -0.07,
                 0.60,  1.00,  0.60, -0.1, -0.07,
                 0.40,  0.60,  1.00, -0.1, -0.07,
                 -0.10, -0.10, -0.10,  1.0,  0.90,
                 -0.07, -0.07, -0.07,  0.9,  1.00), nrow = 5, byrow = T)
  af <- c(0.35, 0.3, 0.4, 0.72, 0.75)

  ## generate genotype data only
  ind1 <- resample_inddata(N = 10,
                           R_LD = list(A),
                           af = af,
                           J = 12)
  expect_equal(nrow(ind1$X), 10)
  expect_equal(ncol(ind1$X), 12)
  expect_identical(ind1$af, c(af, af, af[1:2]))

  # generate some effects
  G <- matrix(0, nrow =2, ncol = 2)
  G[1,2] <- 0.5
  dat <- sim_mv(N = 0,
                J = 12,
                h2 = c(0.03, 0.05),
                pi = 0.5,
                G = G,
                R_LD = list(A),
                af = af)
  N_df <- data.frame("trait_1" = c(1, 1, 0), "trait_2" = c(1, 0, 1),
                     N = c(6, 1, 3))
  ## generate phenotypes only for previous genotype data
  ind2 <- resample_inddata(N = N_df, dat = dat,
                           genos = ind1$X,
                           R_LD = list(A),
                           af = af)
  expect_identical(ind2$Sigma_G, dat$Sigma_G, tolerance = 1e-15)
  expect_identical(ind2$Sigma_E, dat$Sigma_E, tolerance = 1e-15)
  #expect_identical(ind2$X, ind1$X)

  ## generate both phenotypes and genotypes
  ind3 <- resample_inddata(dat = dat,
                           N = N_df,
                           R_LD = list(A),
                           af = af,
                           calc_sumstats = TRUE)
  expect_identical(ind3$Sigma_G, dat$Sigma_G, tolerance = 1e-15)
  expect_identical(ind3$Sigma_E, dat$Sigma_E, tolerance = 1e-15)


  ## change environmental variance
  R_E <- matrix(0.7, nrow = 2, ncol = 2)
  diag(R_E) <- 1
  ind4 <- resample_inddata(dat = dat,
                           N = N_df,
                           R_LD = list(A),
                           af = af,
                           new_env_var = c(0.5, 1.8),
                           new_R_E = R_E)
  expect_equal(diag(ind4$Sigma_E), c(0.5, 1.8))
  ## use data that were on the sd scale
  dat_sd <- rescale_sumstats(dat, output_geno_scale = "sd")
  ind5 <- resample_inddata(dat = dat_sd,
                           N = N_df,
                           R_LD = list(A),
                           af = af)
  expect_equal(ind5$beta_joint, dat$beta_joint)

  ## change LD
  ind6 <- resample_inddata(dat = dat_sd,
                           N = N_df,
                           R_LD = list(matrix(rev(A), nrow = 5)),
                           af = rev(af),
                           new_env_var = c(0.5, 1.8),
                           new_R_E = R_E)
})


