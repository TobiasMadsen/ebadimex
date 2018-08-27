library(testthat)

context("Differential Methylation")

###
# Simulate test datasets
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

# Differential Methylation
grouping <- rep(c(T,F), times=c(100,400))
meth <- simulateModel(param, grouping, expr = expr)

# No differential Methylation
param_null <- param
param_null[8:9] <- 0
meth_null <- simulateModel(param_null, grouping, expr)

###
# Tests
###

test_that("Detect TP 2",{
  log_p <- testDM(meth, expr, grouping)[['p_location']]
  expect_lt(log_p, log(1e-3))
})

test_that("Detect TN 2",{
  log_p <- testDM(meth_null, expr, grouping)[['p_location']]
  expect_gt(log_p, log(5e-2))  
})

test_that("Agree with F-test",{
  # Simulate data
  N1 <- 3
  N2 <- 7
  N <- N1 + N2

  expr <- rnorm(N, 0, 2)
  expr <- expr - mean(expr)
  grouping <- rep(c(T,F), times = c(N1, N2))
  meth <- .1*expr + MASS::mvrnorm(N, mu = c(2,2), Sigma = matrix(c(.18,.06,.06,.22),2,2))
  
  # F-test
  myAnova <- anova(lm(meth ~ expr + grouping))
  f_test <-   myAnova$`Pr(>F)`[3]
  
  # Ebadimex
  ebadimex_test <- testDM(meth = meth, expr = expr, grouping = grouping)
  
  expect_lt(abs(ebadimex_test$p_location - log(f_test)), 1e-9)
})

test_that("Agree with moderated F-test",{
  N1 <- 4
  N2 <- 6
  Lambda <- matrix(c(0.18,0.06,0.06,0.22),2,2)
  nu0 <- 10
  N <- N1 + N2
  
  # Simulate data
  expr <- rnorm(N, 0, 2)
  expr <- expr - mean(expr)
  grouping <- rep(c(T,F), times = c(N1, N2))
  Sigma <- nu0*MCMCpack::riwish(v = nu0, S = Lambda)
  meth <- 1.0*expr + MASS::mvrnorm(N, mu = c(2,2), Sigma = Sigma)
  
  # Moderated F-test by hand
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
  SSD1_tilde <- nu0*Lambda + SSD1
  Q2n_tilde <- det(SSD1_tilde) / det(SSD1_tilde + SSD2)
  
  # Test
  modFTestStat <- (1 - Q2n_tilde) / Q2n_tilde * (N+nu0-3+1-2) / 2
  modFTest <- pf(modFTestStat, 2, nu0+N-3+1-2, lower.tail = F)
 
   
  # Ebadimex
  mp <- list(Lambda = Lambda, nu0 = nu0)
  ebadimex_test <- testDM(meth = meth, expr = expr, grouping = grouping, prior = mp)
  
  expect_lt(abs(ebadimex_test$f_test_statistic - modFTestStat), 1e-9)
  expect_lt(abs(ebadimex_test$p_location - log(modFTest)), 1e-9)
})