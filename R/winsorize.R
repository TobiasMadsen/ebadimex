#' Winsorize Grouped Observations
#' @param x observations to be winsorized
#' @param g logical vector indicating group relationship
#' @param k_win degree of winsorization
#' @param eps convergence tolerance in winsorization procedure
#' @param boost multiply difference from mean by a factor to reduce bias on scale estimate(s)
#' @param pooled use a common standard deviation estimate for each group
#' 
#' @export
winsorize <- function(x, g, k_win, max_iter = 10, eps = 1e-4, boost = FALSE, pooled = TRUE){
  x_star <- x
  n_1 <- sum(g)
  n_2 <- sum(!g)
  
  # Winsorization procedure
  iter <- 0
  while(iter < max_iter){
    iter <- iter + 1

    # Estimate mean in both groups
    m1 <- mean(x_star[g], na.rm = TRUE)
    m2 <- mean(x_star[!g], na.rm = TRUE)
    m <- rep(0, length(x))
    m[g] <- m1
    m[!g] <- m2
    
    # Robust estimate of scale in both groups
    sc <- rep(0, length(x_star))
    if(pooled){
      sc1 <- sc2 <- mean( abs(x_star-m), na.rm = TRUE ) / 0.7738
    } else {
      sc1 <- mean( abs(x_star-m)[g], na.rm = TRUE ) / 0.7738
      sc2 <- mean( abs(x_star-m)[!g], na.rm = TRUE ) / 0.7738
    }
    sc[g] <- sc1
    sc[!g] <- sc2
    x_star > sc * k_win + m + eps
    
    # Correct extreme observations
    has_changed <- any( x_star > sc * k_win + m + eps) | any( x_star < - sc * k_win + m - eps)
    x_star <- pmax( - sc * k_win + m , pmin( sc * k_win + m, x_star ))
    
    # Break if no changes
    if( ! has_changed )
      break
  }

  n_unmodified_1 <- sum(abs(x_star - x)[g] < 1e-9)  
  n_unmodified_2 <- sum(abs(x_star - x)[!g] < 1e-9)

    if(boost){
    # Boost observations. Huber, Robust Statistics, (1.48)
    if( pooled ){
      n_unmodified <- n_unmodified_1 + n_unmodified_2
      n_total   <- n_1 + n_2
      x_star <- m + (n_total / n_unmodified) * (x_star - m)
    } else {
      x_star[g] <- m1 + (n_1 / n_unmodified_1) * (x_star[g] - m1)
      x_star[!g] <- m2 + (n_2 / n_unmodified_2) * (x_star[!g] - m2)
    }
  }
  
  return( list(x_star = x_star, 
               df1 = n_1 - 1 - 2*(n_1 - n_unmodified_1), 
               df2 = n_2 - 1 - 2*(n_2 - n_unmodified_2)) )
}