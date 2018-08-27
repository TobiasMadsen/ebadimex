library(testthat)

context("Winsorization")

test_that("Unchanged observations",{
  x <- c(1:5,1:5)
  g <- rep(c(T,F), each = 5)
  x_star <- winsorize(x, g, k_win = 2)$x_star
  
  expect_equal(x_star, x)
})

test_that("One outlier",{
  x <- c(1:10, 30)
  g <- rep(T, 11)
  x_star <- winsorize(x, g, k_win = 2)$x_star
  
  expect_lt(x_star[11], x[11])
})

test_that("Degree of winsorization",{
  x <- c(1:10, 30)
  g <- rep(T, 11)
  x_star_1 <- winsorize(x, g, k_win = 2)$x_star
  x_star_2 <- winsorize(x, g, k_win = 1.5)$x_star # More agressive winsorization
  
  expect_lt(x_star_2[11], x_star_1[11])
})

test_that("Two groups",{
  x <- c(1:10, 30)
  x <- c(x, x + 5)
  g <- rep(c(T,F), each = 11)
  
  x_star <- winsorize(x, g, k_win = 2)$x_star
  
  expect_lt( abs( x[22] - (x[11] + 5)), 1e-9)
})

test_that("Two groups different variances",{
  x <- c(1:10, 30)
  x <- c(x, 2*x + 5)
  g <- rep(c(T,F), each = 11)
  
  x_star <- winsorize(x, g, k_win = 2, pooled = F)$x_star
  x_star_sub1 <- winsorize(x[1:11], T, k_win = 2, pooled = T)$x_star
  x_star_sub2 <- winsorize(x[12:22], T, k_win = 2, pooled = T)$x_star
  
  expect_lt( sum(abs( x_star[1:11] - x_star_sub1)) , 1e-9)
  expect_lt( sum(abs( x_star[12:22] - x_star_sub2)) , 1e-9)
})


test_that("Outlier border",{
  x <- c(rep(0, 10), 40)
  g <- rep(T, 11)
  
  # Outlier border
  border <- 11/2 * 0.7738
  
  x_star_1 <- winsorize(x, g, k_win = border + 1e-1)$x_star # Inside border
  x_star_2 <- winsorize(x, g, k_win = border - 1e-1)$x_star # Outside border
  
  expect_equal(x_star_1[11], 40)
  expect_lt(x_star_2[11], 40)
})
