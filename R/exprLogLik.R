#' Compute Robust-Normal likelihood
#' 
robustNormLogLik <- function(x, m, sd, k_win){
  - log(2*pi) / 2 - log(sd) - ifelse( abs(x - m) < k_win * sd, (x - m)^2, 2 * k_win * sd * abs(x - m) - (k_win*sd)^2 ) / 2 / sd^2
}

#' Calculate Log-Likelihood-Ratio for Expression data
#' 
#' To ensure robustness the log-lik is linearized when
#' the expression measurement is more than k_win standard
#' deviations away from the group mean.
#' 
#' @export
exprLogLik <- function(expr, param, same_sd = T, k_win = 2.5){
  if(same_sd){
    loglik1 <- robustNormLogLik(expr, param$m1, param$sd_reg, k_win)
    loglik2 <- robustNormLogLik(expr, param$m2, param$sd_reg, k_win)
  } else {
    loglik1 <- robustNormLogLik(expr, param$m1, param$sd1_reg, k_win)
    loglik2 <- robustNormLogLik(expr, param$m2, param$sd2_reg, k_win)
  }
  # Log-likelihood ratio
  loglik1 - loglik2
}