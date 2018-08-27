library(testthat)

context("Kuhlback Leibler Divergence")

test_that("Univariate",{
  expect_lt( abs(kl_normal(0.1, 2, 0.1, 2)) , 1e-16)
})

test_that("Bivariate",{
  mu1 <- c(1,2)
  Sigma1 <- matrix(c(4,1,1,4),2,2)

  expect_lt( abs(kl_multivariate(mu1, Sigma1, mu1, Sigma1)), 1e-16)
})