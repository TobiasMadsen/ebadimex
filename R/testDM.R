
#' Test differential Methylation
#' 
#' @export
testDM <- function(meth, expr, grouping, win = F, prior = NULL){
  fm <- fitMeth(meth, expr, grouping, win = win, prior = prior)

  N <- nrow(meth)
  n1 <- sum(grouping)
  n2 <- sum(!grouping)
  
  # Moderated Bartlett test
  f1 <- n1 - 2 + fm$nu0
  f2 <- n2 - 2
  SSDgT <- fm$SSDg1 + fm$nu0 * fm$Lambda_0
  SSDgF <- fm$SSDg2
  p_obs_bartlett <- bartlettCov(SSD1 = SSDgT, SSD2 = SSDgF, f1 = f1, f2 = f2)
  
  # Moderated f-statistic    
  SSD1 <- fm$SSD1
  SSD2 <- fm$SSD2
  Q2n <- det(SSD1 + fm$nu0 * fm$Lambda_0) / det(fm$nu0 * fm$Lambda_0 + SSD1 + SSD2)
  nu0 <- fm$nu0

  f_test_statistic <- (1 - Q2n) / Q2n * (N+nu0-3+1-2) / 2
  p_obs_location <- pf(f_test_statistic, df1 = 2, df2 = N+nu0-3+1-2, lower.tail = F, log.p = T)
  
  # Compute Kuhlback-Leibler Divergences
  kl_12 <- kl_multivariate(c(0,0), fm$param$Sigma_g1_reg, c(fm$param$diff_x, fm$param$diff_y), fm$param$Sigma_g2_reg)
  kl_21 <- kl_multivariate(c(0,0), fm$param$Sigma_g2_reg, c(fm$param$diff_x, fm$param$diff_y), fm$param$Sigma_g1_reg)

  return(list(p_location = p_obs_location,
              f_test_statistic = f_test_statistic,
              p_bartlett = p_obs_bartlett,
              kl_12 = kl_12,
              kl_21 = kl_21))
}