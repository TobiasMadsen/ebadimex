#' Calculate Log-Likelihood-Ratio for Expression data
#' 
#' @export
methLogLik <- function(meth, expr, param, same_cov = T){
  # Determine probe status: Promoter or Genebody
  N <- nrow(meth)

  # Calculate likelihood of each sample
  M1 <- as.matrix(cbind(expr, 1, 1)) %*% param$alpha
  M2 <- as.matrix(cbind(expr, 1, 0)) %*% param$alpha

  if(same_cov){  
    loglik1 <- mvtnorm::dmvnorm(meth-M1, mean = c(0,0), sigma = param$Sigma_reg)
    loglik2 <- mvtnorm::dmvnorm(meth-M2, mean = c(0,0), sigma = param$Sigma_reg)
  } else {
    loglik1 <- mvtnorm::dmvnorm(meth-M1, mean = c(0,0), sigma = param$Sigma_g1_reg)
    loglik2 <- mvtnorm::dmvnorm(meth-M2, mean = c(0,0), sigma = param$Sigma_g2_reg)
  }
  # Return loglik ratio
  loglik1 - loglik2
}
