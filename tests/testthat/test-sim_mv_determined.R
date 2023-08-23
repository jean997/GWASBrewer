test_that("sim_mv_determined works", {
  G <- matrix(c(0, 0.5, 0, 0), nrow = 2, byrow =T)
  my_effects <- matrix(0, nrow = 10, ncol = 2)
  my_effects[c(1, 5),1] <- c(-0.008, 0.01)
  my_effects[c(3, 6, 9), 2] <- c(-0.02, 0.06, 0.009)

  set.seed(11111)
  sim_dat1 <- sim_mv_determined(N = 0,
                                 direct_SNP_effects_joint = my_effects,
                                 geno_scale = "sd",
                                 pheno_sd = 1,
                                 G=G,
                                 est_s = TRUE)
  expect_equal(sim_dat1$direct_SNP_effects_joint, my_effects)


  # now give effects as per-allele
  af <- rbeta(10, 1, 5)
  sim_dat2 <- sim_mv_determined(N = 0,
                               direct_SNP_effects_joint = my_effects,
                               geno_scale = "allele",
                               pheno_sd = 1,
                               G=G,
                               est_s = TRUE,
                               af = af)
  expect_equal(sim_dat1$direct_SNP_effects_joint, my_effects)

  # give effects as per-sd but proved af so converted to per-allele
  sim_dat3 <- sim_mv_determined(N = 0,
                                direct_SNP_effects_joint = my_effects,
                                geno_scale = "sd",
                                pheno_sd = 1,
                                G=G,
                                est_s = TRUE,
                                af = af)
  expect_equal(sim_dat3$direct_SNP_effects_joint*sqrt(2*af*(1-af)), my_effects)

  # change pheno sd
  sim_dat4 <- sim_mv_determined(N = 0,
                                direct_SNP_effects_joint = my_effects,
                                geno_scale = "sd",
                                pheno_sd = c(0.5, 0.8),
                                G=G,
                                est_s = TRUE,
                                af = af)
  expect_equal(sim_dat4$direct_trait_effects[1,2], 0.5*0.5/0.8)
  A <- kronecker(sqrt(2*af*(1-af)), c(0.5, 0.8)) |> matrix(nrow = 10, byrow = T)
  expect_equal(sim_dat4$direct_SNP_effects_joint*A, my_effects)
})
