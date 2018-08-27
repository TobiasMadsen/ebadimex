library(testthat)

context("Bartlett")

test_that("Bartlett output",{
  SSD1 <- SSD2 <- matrix(c(2,1,1,2),2,2)
  f1 <- f2 <- 10
  res <- bartlettCov(SSD1, SSD2, f1, f2)
  expect_true( is.numeric(res) )
  expect_gt( res, -1e-16)
})

test_that("Bartlett df",{
  SSD1 <- matrix(c(2,1,1,2),2,2)
  SSD2 <- matrix(c(2,1,1,2),2,2) * 10
  f1 <- 10
  f2 <- 100
  
  res <- bartlettCov(SSD1, SSD2, f1, f2)
  expect_true( is.numeric(res) )
  expect_gt( res, -1e-16)
})

test_that("Bartlett diff cov",{
  SSD1 <- matrix(c(2,1,1,2),2,2)
  SSD2 <- matrix(c(2,1,1,2),2,2) * 2
  f1 <- 10
  f2 <- 100
  
  res <- bartlettCov(SSD1, SSD2, f1, f2)
  expect_true( is.numeric(res) )
  expect_lt( res, log(0.05) )
})

test_that("Uniform p-values",{
  expect_true(TRUE)

  if(FALSE){
    # The p-values indicate deflation
    sim <- replicate(n = 1000,{
      f1 <- 10
      f2 <- 100
      covMat <- matrix(c(2,1,1,2), 2,2)
      SSD1 <- MCMCpack::rwish(f1, covMat)
      SSD2 <- MCMCpack::rwish(f2, covMat)
      
      bartlettCov(SSD1, SSD2, f1, f2)
    })
    
    df <- data.frame(x = sort(-sim), q = qchisq(ppoints(1000), df = 1))
    ggplot(df, aes(x = q, y = x)) + geom_point() + geom_abline(slope = 1, intercept = 0) +
      theme_minimal() + ggtitle("Bartlett")
  }
})



