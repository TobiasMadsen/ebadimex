library(testthat)
library(numDeriv)

context("Expr Log Lik")

test_that("Agree with dnorm",{
  # Agrees with dnorm in normal neighborhood
  x <- 0.25
  m <- 0.15
  sd <- 0.3
  expect_lt( abs( robustNormLogLik(x, m, sd, 2.5) - dnorm(x, m, sd, log = T) ), 1e-9)
  
  x <- -0.25
  m <- 0.15
  sd <- 0.4
  expect_lt( abs( robustNormLogLik(x, m, sd, 2.5) - dnorm(x, m, sd, log = T) ), 1e-9)
})

test_that("Agree with dnorm two groups",{
  x <- rnorm(10, 0, 1)
  y <- rnorm(10, 1, 0.5)
  irksome_loglik <- exprLogLik(c(x,y), param = list(m1 = 0, m2 = 1, sd_reg = 0.75), k_win = 10)
  dnorm_loglik <- dnorm(c(x,y), 0, 0.75, log = T) - dnorm(c(x,y), 1, 0.75, log = T)
  
  expect_lt(max(abs(irksome_loglik - dnorm_loglik)), 1e-9)
})

test_that("Agree with dnorm two groups different variances",{
  x <- rnorm(10, 0, 1)
  y <- rnorm(10, 1, 0.5)
  irksome_loglik <- exprLogLik(c(x,y), param = list(m1 = 0, m2 = 1, sd_reg = 0.75, sd1_reg = 0.5, sd2_reg = 1), same_sd = F, k_win = 10)
  dnorm_loglik <- dnorm(c(x,y), 0, 0.5, log = TRUE) - dnorm(c(x,y), 1, 1, log = TRUE)
  
  expect_lt(max(abs(irksome_loglik - dnorm_loglik)), 1e-9)
})

test_that("Continuous",{
  # Likelihood function is continuous
  expect_true( all( diff(robustNormLogLik(seq(0, 5, 0.01), 0, 0.5, 2.5)) < 0 ))
})

test_that("Linear outside normal neighborhood",{
  diff_lik <- diff(robustNormLogLik(seq(1,3,0.01), 0, 0.5, 2))
  expect_lt( max(abs(diff_lik - diff_lik[1])) , 1e-9)
})

test_that("Sign of loglik-ratio",{
  # loglik_true - loglik_false
  x <- rnorm(100, 0, 1)
  y <- rnorm(100, 10, 1)
  expr <- c(x,y)
  grouping <- rep(c(T,F), times = c(100, 100))
  
  fm <- fitExpr(expr, grouping, k_win = 2.5)
  ell <- exprLogLik(expr, fm, k_win = 2.5)
  
  expect_true( min(ell[1:100]) > max(ell[101:200]))
})
