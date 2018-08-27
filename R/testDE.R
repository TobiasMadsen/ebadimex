##################################################
# Fit expression and evaluate differential expression
##################################################

#' Test for differential Expression
#' 
#' Fit expression model and evaluate differential expression
#' Robust regression model with prior on standard deviation
#' @param grouping Logical vector indicating group relationship 
#' @param expr   Vector with expression values
#' @param prior  
#' @export
testDE <- function(expr, grouping, prior = NULL, k_win = 2.5, pooled = F, winz_df = T){
  fe <- fitExpr(expr, grouping, prior, k_win, pooled = pooled)
  n  <- length(grouping)
  n1 <- sum(grouping)  
  n2 <- sum(!grouping)
  
  # Asses variance difference
  # Use regularized estimate for smaller group
  df1 <- ifelse(winz_df, fe$df1, n1 - 1)
  df2 <- ifelse(winz_df, fe$df2, n2 - 1)
  if( sum(grouping) > sum(!grouping)){ # TRUE Group is larger
    f <- (fe$sd1)^2 / (fe$sd2_reg)^2
    df2 <- df2 + fe$d0
  } else {
    f <- (fe$sd1_reg)^2 / (fe$sd2)^2
    df1 <- df1 + fe$d0
  }

  p_obs_f <- 
    log(2) + ifelse(f > 1, 
                    pf(f, df1 = df1, df2 = df2, lower.tail = F, log.p = T), 
                    pf(1/f, df1 = df2, df2 = df1, lower.tail = F, log.p = T))
  
  # Asses Differential Expression using moderated t-statistic
  t <- (fe$m1 - fe$m2) / fe$sd_reg / sqrt(1/sum(grouping)+1/sum(!grouping))
  p_obs_t <- log(2) + pt(abs(t), df = fe$d0+length(grouping)-2, log.p = TRUE, lower.tail = FALSE)
  
  # Asses Differential Expression using moderated t-statistic unequal variances
  # Moderated Welch T-test
  t_welch <- (fe$m1 - fe$m2) / sqrt(fe$sd1_reg^2/n1 + fe$sd2_reg^2/n2)
  df_welch <- (fe$sd1_reg^2/n1 + fe$sd2_reg^2/n2)^2 / (fe$sd1_reg^4 / n1^2 / (n1-1+fe$d0) + fe$sd2_reg^4 / n2^2 / (n2-1+fe$d0))
  p_obs_t_welch <- log(2) + pt(abs(t_welch), df = df_welch, log.p = TRUE, lower.tail = FALSE)
  
  # Compute Kuhlback-Leibler Divergences
  kl_12 <- kl_normal(fe$m1, fe$sd1_reg, fe$m2, fe$sd2_reg)
  kl_21 <- kl_normal(fe$m2, fe$sd2_reg, fe$m1, fe$sd1_reg)
  
  # Return
  return( list(p_location = p_obs_t, 
               p_location_welch = p_obs_t_welch, 
               p_scale = p_obs_f,
               kl_12 = kl_12,
               kl_21 = kl_21
  ))
}
