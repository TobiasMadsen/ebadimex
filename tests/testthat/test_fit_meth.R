library(testthat)

context("Fit Methylation")

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
           0.3,
           0.2)

grouping <- rep(c(T,F), times=c(100,400))
meth <- simulateModel(param, grouping, expr)

###
# Tests
###

test_that("Identifies parameters fitMeth",{
  fm <- fitMeth(meth, expr, grouping)
  
  # Slopes
  expect_lt( abs(fm$param$a_pr - (-1)), 1e-2)
  expect_lt( abs(fm$param$a_gb - (-0.5)), 1e-2)
  
  # Difference between groups
  expect_lt( abs(fm$param$diff_x - 0.3), 1e-1)
  expect_lt( abs(fm$param$diff_y - 0.2), 1e-1)

  # Returns SSD Matrix
  expect_false( is.null(fm$SSD1) )
  
  # Covariance matrix
  expect_lt( sum(abs( fm$SSD1/500 - matrix(c(0.2+exp(-1.6)/48,0.1,0.1,0.3+exp(-1.8)/10),2,2)) ), 2e-1)
  
  # Returns loglik
  expect_type( fm$loglik, "double")
})

test_that("Estimation of SSD1 and SSD2",{
  # Simulate data
  N1 <- 3
  N2 <- 7
  N <- N1 + N2
  
  expr <- rnorm(N, 0, 2)
  expr <- expr - mean(expr)
  grouping <- rep(c(T,F), times = c(N1, N2))
  meth <- .1*expr + MASS::mvrnorm(N, mu = c(2,2), Sigma = matrix(c(.18,.06,.06,.22),2,2))
  
  # Design matrix
  Z1 <- matrix( c(expr, 
                  rep(c(1,0), times = c(N1, N2)), 
                  rep(c(0,1), times = c(N1, N2))),
                N, 3)
  Z2 <- matrix( c(expr,
                  rep(1, times = N)),
                N, 2)
  
  # SSD Matrices
  SSD1 <- t(meth) %*% meth - t(meth) %*% Z1  %*% solve(t(Z1) %*% Z1 ) %*% t(Z1) %*% meth
  SSD2 <- (t(meth) %*% Z1  %*% solve(t(Z1) %*% Z1 ) %*% t(Z1) %*% meth) -  
    (t(meth) %*% Z2  %*% solve(t(Z2) %*% Z2 ) %*% t(Z2) %*% meth)
  
  # Ebadimex fit
  fm <- fitMeth(meth, expr, grouping)
  fm$SSD1
  fm$SSD2
  
  expect_lt(sum(abs(fm$SSD1 - SSD1)), 1e-9)
  expect_lt(sum(abs(fm$SSD2 - SSD2)), 1e-9)
})

# Covariance Matrix is Symmetric
test_that("Covariance Matrix Symmetric",{
  fm <- fitMeth(meth, expr, grouping)
  
  expect_true( isSymmetric(fm$param$Sigma) )
  expect_true( isSymmetric(fm$SSD1) )
})