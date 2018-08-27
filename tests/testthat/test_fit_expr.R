library(testthat)

context("Fit Expression")

test_that("Returns list of parameters",{
  x <- rnorm(30, 0, 0.4)
  y <- rnorm(50, 0.2, 0.4)
  
  fe <- fitExpr(c(x,y), grouping = rep(c(T,F), times = c(30, 50)), k_win = Inf)
  
  expect_true( is.list(fe) )
  expect_true( all(c("m1","m2","sd","sd_reg","sd1_reg","sd2_reg","d0","s0") %in% names(fe)) )
})

test_that("Agree with lm without winsorization",{
  x <- rnorm(30, 0, 0.4)
  y <- rnorm(50, 0.2, 0.4)
  
  fe <- fitExpr(c(x,y), grouping = rep(c(T,F), times = c(30, 50)), k_win = Inf)
  
  df <- data.frame(value = c(x,y), group = rep(c(T,F), times = c(30, 50)))
  mod <- lm(value ~ group + 0, data = df)
  
  expect_lt( abs(coef(mod)['groupTRUE'] - fe$m1), 1e-9)  
  expect_lt( abs(coef(mod)['groupFALSE'] - fe$m2), 1e-9)  
  
  expect_lt( abs(sigma(mod) - fe$sd), 1e-9)
})

