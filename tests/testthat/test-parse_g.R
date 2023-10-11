# Define the tests for the G_from_text function
test_that("G_from_text parses text input correctly", {
  # Define some example text input
  text <- "X <- Y 0.25
           A <- X 0.1
           A <- Y sqrt(0.34)"
  
  # Call the G_from_text function
  G_matrix <- G_from_text(text)
  
  # Define the expected output G matrix (you may need to adjust this based on your function's logic)
  expected_G_matrix <- matrix(c(0, 0.25, sqrt(0.34), 0, 0, 0.1, 0, 0, 0), nrow = 3, ncol = 3,
                              dimnames = list(c("Y", "X", "A"), c("Y", "X", "A")), byrow=T)
  
  # Check if the generated G matrix matches the expected G matrix
  expect_equal(G_matrix, expected_G_matrix)
})

test_that("G_from_text handles invalid input gracefully", {
  # Define some invalid input
  invalid_text <- "X -> Y 0.25\nInvalidLine\nA <- Z 0.1"
  
  # Check if the function throws an error when given invalid input
  expect_error(G_from_text(invalid_text), "Expect there to be four tokens, space separated")
})



test_that("G_from_df generates G matrix from a valid data frame", {
  # Define a sample data frame
  df <- data.frame(
    i = c("Y", "X", "X", "A", "Y", "A"),
    j = c("X", "A", "Y", "X", "A", "Y"),
    eff = c(0.25, 0.24, 0.34, 0.2, 0.15, 0.3)
  )
  
  # Call the G_from_df function
  G_matrix <- G_from_df(df)
  
  # Define the expected output G matrix (you may need to adjust this based on your function's logic)
  expected_G_matrix <- matrix(c(0, 0.25, 0.15, 0.34, 0, 0.24, 0.3, 0.2, 0), nrow = 3, ncol = 3,
                              dimnames = list(c("Y", "X", "A"), c("Y", "X", "A")), byrow=T)
  
  # Check if the generated G matrix matches the expected G matrix
  expect_equal(G_matrix, expected_G_matrix)
})

test_that("G_from_df handles invalid input gracefully", {
  # Define an empty data frame
  empty_df <- data.frame()
  
  # Check if the function throws an error when given an empty data frame
  expect_error(G_from_df(empty_df), "nrow")
  
  # Define a data frame with missing columns
  missing_cols_df <- data.frame(i = c("Y", "X"), eff = c(0.25, 0.24))
  
  # Check if the function throws an error when given a data frame with missing columns
  expect_error(G_from_df(missing_cols_df), "all")
})


