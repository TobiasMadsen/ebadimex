
logDiff <- function(a,b){
  # Computes log(exp(a) - exp(b))
  # Numerically stable
  
  if(is.infinite(b) && sign(b) == -1 )
    return(a)

  log(exp(a-max(a,b)) - exp(b-max(a,b))) + max(a,b)
}

logSum <- function(a,b){
  # Computes log(exp(a)+exp(b))
  
  if(is.infinite(a) && sign(a) == -1 )
    return(b)
  if(is.infinite(b) && sign(b) == -1 )
    return(a)
  
  log(exp(a-max(a,b)) + exp(b-max(a,b))) + max(a,b)
}

#' Bartlett Test of Covariance Homogeneity
#' @export
bartlettCov <- function(SSD1, SSD2, f1, f2){
  logQ1 <- f1/2*log(abs(det(SSD1))) + f2/2*log(abs(det(SSD2))) - (f1+f2)/2*log(abs(det(SSD1+SSD2))) +
    (f1+f2)*log(f1+f2) - f1*log(f1) - f2*log(f2)
  
  rho <- 1 - (1/f1+1/f2-1/(f1+f2))*(13/18) # (2*p^2+3*p-1)/(6*(p+1)*(q-1))
    
  z <- -2*rho*logQ1
  
  omega2 <- (6*(4*(1/f1^2+1/f2^2-1/(f1+f2)^2)-6*(1-rho)^2 )) / 48 / rho^2
  k <- 3 # Reduction in free parameters

  # p_obs_bartlett <- log(pchisq(z, k, lower.tail = F) - omega2*(pchisq(z, df = k, lower.tail = F) - pchisq(z, df = k + 4, lower.tail = F)))
  # p_obs_bartlett <- log1p( -(pchisq(z, df = k)+omega2*(pchisq(z, df = k + 4) - pchisq(z, df = k))) )
  # Numerically stable calculation of the above
  tmp1 <- pchisq(z, k, lower.tail = FALSE, log.p = TRUE)
  tmp2 <- logDiff(pchisq(z, df = k + 4, lower.tail = FALSE, log.p = TRUE), pchisq(z, df = k, lower.tail = FALSE, log.p = TRUE)) 
  return( ifelse(omega2 > 0, logSum(tmp1, tmp2 + log(omega2)), logDiff(tmp1, tmp2 + log(-omega2))) )

}
