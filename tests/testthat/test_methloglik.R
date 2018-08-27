library(testthat)

context("Meth Log Lik")

###
# Simulated test dataset
###

set.seed(1)
expr <- rnorm(500, -9, 2)
param <- c(-1,
           -0.5,
           -2.6,
           -0,
           covToPar(matrix(c(0.2,0.1,0.1,0.3),2,2)),
           -1.6,
           -1.8,
           1.7,
           1.5)

grouping <- rep(c(T,F), times=c(100,400))
meth <- simulateModel(param, grouping, expr = expr)

fm <- fitMeth(meth, expr, grouping)

# Takes output from fitmeth
test_that("Takes output from fitMeth",{
  # Full methylation
  loglik <- methLogLik(meth, expr, fm$param)
  
  expect_equal(length(loglik), 500)
  expect_true( ! any(is.na(loglik)))
  expect_true( ! any(is.infinite(loglik)))
})

test_that("Takes output from fitMeth Different Covariances",{
  # Full methylation
  loglik <- methLogLik(meth, expr, fm$param, same_cov = T)
  
  expect_equal(length(loglik), 500)
  expect_true( ! any(is.na(loglik)))
  expect_true( ! any(is.infinite(loglik)))
})

# Has correct sign
test_that("Sign of loglik-ratio",{
  # loglik_true - loglik_false
  loglik <- methLogLik(meth, expr, fm$param)
  expect_gt( quantile(loglik[1:100], probs = 0.1), quantile(loglik[101:500], probs = 0.9))
})
