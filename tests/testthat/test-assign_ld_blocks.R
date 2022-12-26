test_that("assign_ld_blocks works", {
  # J is less than total
  t1 <- assign_ld_blocks(c(2, 4, 5), J = 10)
  expect_length(t1$index, 10)
  expect_equal(t1$last_block_info, c(3, 4))
  expect_equal(t1$rep, rep(1, 3))
  expect_equal(t1$block_index, c(1, 2,4))

  # J is equal to total
  t2 <- assign_ld_blocks(c(2, 4, 5), J = 11)
  expect_length(t2$index, 11)
  expect_null(t2$last_block_info)
  expect_equal(t2$rep, rep(1, 3))
  expect_equal(t2$block_index, c(1, 2,3))

  # J greater than total
  t3 <- assign_ld_blocks(c(2, 4, 5), J = 20)
  expect_length(t3$index, 20)
  expect_equal(t3$last_block_info, c(3, 3))
  expect_equal(t3$rep, c(rep(1, 3), rep(2, 3)))
  expect_equal(t3$block_index, c(1, 2,3,1,2,4))
})
