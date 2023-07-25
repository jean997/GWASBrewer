test_that("sim_mv", {
    G <- matrix(c(0, sqrt(0.25), 0, sqrt(0.15), 
                0, 0, 0, sqrt(0.1), 
                sqrt(0.2), 0, 0, -sqrt(0.3), 
                0, 0, 0, 0), nrow = 4, byrow = TRUE)
    colnames(G) <- row.names(G) <- c("X", "Y", "Z", "W")

    x <- sim_mv(G = G,
                    J = 50000,
                    N = 60000,  
                    h2 = c(0.3, 0.3, 0.5, 0.4), 
                    pi = 1000/50000) %>% suppressWarnings()

    expect_true(all(!is.na(x$beta_hat)))
    expect_true(sum(x$beta_hat == 0) / length(x$beta_hat) < 0.01)
})
