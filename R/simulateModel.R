##################################################
# Simulate a dataset
##################################################

#' @export
simulateModel <- function(par, group, expr){
  # Parameters
  ax      <- par[1]     
  ay      <- par[2]     
  bx      <- par[3]     
  by      <- par[4]     
  th_1    <- par[5]
  th_2    <- par[6]  
  th_3    <- par[7]
  diff_x  <- par[8] 
  diff_y  <- par[9]
  
  sim <- mvtnorm::rmvnorm(length(group), mean = c(0,0), sigma = parToCov(c(th_1,th_2,th_3)))
  sim <- sim + c(expr*ax+bx, expr*ay+by) + c(group * diff_x, group * diff_y)
  
  return(sim)
}