test_that("fast_lm is equivalent to lm", {
  X <- matrix(rbinom(n = 100*2, size = 2, prob = 0.3), nrow  = 100)
  Y <- matrix(rnorm(n = 100), nrow = 100)
  fast_lm_res <- fast_lm(X, Y)
  lm_res <- purrr::map_dfr(1:2, function(i){
    f <- lm(Y ~ X[,i])
    data.frame("bhat_1" = f$coefficients[2], "s_1" = sqrt(vcov(f)[2,2]) )
  })
  expect_equal(lm_res$bhat_1, fast_lm_res$bhat_1)
  expect_equal(lm_res$s_1, fast_lm_res$s_1)
})
