library(testthat)

context("Differential Expression")

# Agree with t-test
test_that("Agreement with t.test",{
  x <- rnorm(30, 0, 0.4)
  y <- rnorm(50, 0.2, 0.4)
  
  irksome_test <- testDE(c(x,y), rep(c(T,F), times = c(30, 50)), k_win = 100)[1]
  
  t_test <- log(t.test(x, y, var.equal = T)$p.value)

  
  expect_lt(abs(irksome_test[['p_location']] - t_test), 1e-9)
})

test_that("Agreement with moderated t-test",{
  x <- rnorm(3, 0)
  y <- rnorm(7, 0)
  d0 <- 4
  s0 <- 1.2  
  ep <- cbind(expr = c(1,2), d0 = c(d0,d0), s0 = c(s0,s0))
  
  # Moderated t-test within ebadimex
  ebadimex_test <- testDE(c(x,y), rep(c(T,F), times = c(3, 7)), prior = ep, k_win = 100)[1]
  
  # Compute moderated t by hand
  ssd1 <- sum((x - mean(x))^2)
  ssd2 <- sum((y - mean(y))^2)
  d_g <- 3 + 7 - 2
  s2 <- (ssd1+ssd2) / d_g
  s2_mod <- (d_g * s2 + d0 * s0) / (d_g + d0)
  mod_t_test_statistic <- (mean(x) - mean(y)) / sqrt(s2_mod) / sqrt(1/3+1/7)
  mod_t_test <- log(2*pt(-abs(mod_t_test_statistic), d_g + d0))

  expect_lt(abs(ebadimex_test[['p_location']] - mod_t_test), 1e-9)
})

test_that("Agreement with Welch t.test",{
  x <- rnorm(30, 0, 0.4)
  y <- rnorm(50, 0.2, 0.4)
  
  irksome_test <- testDE(c(x,y), rep(c(T,F), times = c(30, 50)), k_win = 100)
  
  t_test <- log(t.test(x, y, var.equal = F)$p.value)  
  
  expect_lt(abs(irksome_test[['p_location_welch']] - t_test), 1e-9)
})

# Agree with F-test
test_that("Agreement with var.test",{
  x <- rnorm(30, 0, 0.4)
  y <- rnorm(50, 0.2, 0.4)
  
  irksome_test <- testDE(c(x,y), rep(c(T,F), times = c(30, 50)), k_win = 100)
  
  var_test <- log(var.test(x,y)$p.value)
  
  expect_lt(abs(irksome_test[['p_scale']] - var_test), 1e-9)
})

# Agree with moderated F-test
test_that("Agreement with moderated F-test",{
  
})