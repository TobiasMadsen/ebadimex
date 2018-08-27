##################################################
# Unconstrained parameterization 
##################################################

#' @export
parToCov <- function(th){
  th_1 <- th[1]
  th_2 <- th[2]
  th_3 <- th[3]
  
  l_11 <- exp(th_1)
  l_21 <- exp(th_2)
  l_22 <- pi*exp(th_3) / (1+exp(th_3))
  
  Sigma11 <- l_11*l_11;
  Sigma22 <- l_21*l_21;
  Sigma12 <- l_11*l_21 * cos(l_22);
  
  matrix(c(Sigma11,Sigma12,Sigma12,Sigma22), 2,2)
}

#' @export
covToPar <- function(sigma){
  l_11 <- sqrt(sigma[1,1])
  l_21 <- sqrt(sigma[2,2])
  l_22 <- acos(sigma[1,2]/l_11/l_21)
  
  th_1 <- log(l_11)
  th_2 <- log(l_21)
  th_3 <- log(l_22 / (pi - l_22))
  
  return(c(th_1, th_2, th_3))
}
