context("Prior")

##################################################
# Generation of data for tests
##################################################

if(! file.exists("../testdata/data_prior.RData")){
  d <- 4
  nu0 <- 6
  Lambda <- matrix(c(5,1,1,2), 2,2)
  
  covMats <- lapply(1:20000, FUN = function(i){
    Sigma <- MCMCpack::riwish(nu0, Lambda*nu0)
    MCMCpack::rwish(d, Sigma) / d
  })
  save(covMats, file = "ebadimex/tests/testdata/data_prior.RData")
}
load("../testdata/data_prior.RData")

##################################################
# Tests
##################################################

test_that("Infer inverse chi square 1",{
  d <- 4
  nu0 <- 6
  s0 <- 5
  s2 <- s0*rf(20000, d, nu0)
  
  fit <- fitInverseChiSquare(s2, d)  
  
  expect_lt(abs(fit$s0 - s0), 0.2)
  expect_lt(abs(fit$nu0 - nu0), 0.2 )
})

test_that("Infer inverse chi square 2",{
  d <- 4
  nu0 <- 6
  s0 <- 5
  sigma <- s0*nu0/rchisq(20000,nu0)
  s2 <- sigma*rchisq(20000, df = d)/d
  
  fit <- fitInverseChiSquare(s2, d)
  
  expect_lt(abs(fit$s0 - s0), 0.2)
  expect_lt(abs(fit$nu0 - nu0), 0.3 )
})

test_that("Infer inverse Wishart",{
  d <- 4
  nu0 <- 6
  s0 <- matrix(c(5,1,1,2), 2,2)
  
  #! Note covMats generated in beginning of file
  s2 <- t(sapply(covMats, diag))
  
  # Marginal inference
  fit <- fitInverseChiSquare(s2[,1], d)
  expect_lt(abs(fit$nu0 - (nu0-1)), 0.2)
  expect_lt(abs(fit$s0 - (s0[1,1] * nu0 / (nu0-1)) ), 0.2)
  
  fit <- fitInverseChiSquare(s2[,2], d)
  expect_lt(abs(fit$nu0 - (nu0-1)), 0.2)
  expect_lt(abs(fit$s0 - (s0[2,2] * nu0 / (nu0-1)) ), 0.2)
  
  # Joint Inference
  fit <- fitInverseChiSquare(s2, d)
  expect_lt( abs(fit$nu0 - nu0), 0.2)
  expect_lt( abs(fit$s0[1] - s0[1,1]), 0.2)
  expect_lt( abs(fit$s0[2] - s0[2,2]), 0.2)
})

test_that("Fit Prior",{
  d <- 4
  nu0 <- 6
  Lambda <- matrix(c(5,1,1,2), 2,2)
 
  #! covMats generated in beginning of file
  fit <- fitInverseWishart(covMats, d)
  
  expect_lt( abs(fit$nu0 - nu0), 0.2)
  expect_lt( sum(abs(fit$Lambda - Lambda)), 0.1 * 4)
})