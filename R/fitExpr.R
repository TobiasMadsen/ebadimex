
#' Fit expression data
#' @param expr vector with expression values
#' @param grouping logical vector indicating group relationship
#' @param k_win amount of winsorization cap values more than k_win standard deviations away from mean
#' @param pooled use a common standard deviation estimate for each group
#' @export
fitExpr <- function(expr, grouping, prior = NULL, k_win = 2.5, pooled = T){
  K <- length(expr)
  
  # Winsorize expression
  winz_expr <- winsorize(expr, grouping, k_win, boost = T, pooled = pooled)
  expr_star <- winz_expr$x_star
  
  # Compute mean in each group
  m1 <- mean(expr_star[grouping])
  m2 <- mean(expr_star[!grouping])
  
  m <- rep(0, K)
  m[grouping] <- m1
  m[!grouping] <- m2
  
  # Compute standard deviation
  sd <- sqrt( sum((expr_star - m)^2) / (K-2) )
  sd1 <- sqrt( sum( ((expr_star - m)^2)[grouping] ) / (sum(grouping)-1) )
  sd2 <- sqrt( sum( ((expr_star - m)^2)[!grouping] ) / (sum(!grouping)-1) )

  # Apply Prior 
  s0 <- 0
  d0 <- 0
  if(!is.null(prior)){
    m12 <- mean(expr_star)
    
    s0 <- approx(prior[,"expr"], prior[,"s0"], xout = m12, rule = 2)$y
    d0 <- approx(prior[,"expr"], prior[,"d0"], xout = m12, rule = 2)$y    
  }
  
  sd_reg  <- sqrt( (d0*s0 + (K-2)*sd^2) / (d0 + K-2) )
  sd1_reg <- sqrt( (d0*s0 + (sum(grouping)-1)*sd1^2) / (d0 + sum(grouping)-1) )
  sd2_reg <- sqrt( (d0*s0 + (sum(!grouping)-1)*sd2^2) / (d0 + sum(!grouping)-1) )
  
  # Parameters:
  #  - Upper and lower cut values for winsorization
  #  - Means and standard deviation(s)
  list(m1 = m1, # Mean for grouping: TRUE
       m2 = m2, # Mean for grouping: FALSE
       sd = sd, # Standard deviation
       sd_reg = sd_reg,
       sd1 = sd1,
       sd1_reg = sd1_reg,
       sd2 = sd2,
       sd2_reg = sd2_reg,
       d0 = d0,
       s0 = s0,
       df1 = winz_expr$df1,
       df2 = winz_expr$df2) 
}
