#' KL Univariate Normal
#' Compute KL-divergence between two
#' univariate normal distributions
#' @export
kl_normal <- function(mu1, sigma1, mu2, sigma2){
  log(sigma2 / sigma1) + (sigma1^2-sigma2^2+(mu1-mu2)^2)/sigma2^2/2
}

#' KL Multivariate Normal
#' Compute KL-divergence between two
#' multivariate normal distributions
#' @export
kl_multivariate <- function(mu1, Sigma1, mu2, Sigma2){
  c((log(abs(det(Sigma2)/det(Sigma1))) - length(mu1) + sum(diag(solve(Sigma2) %*% Sigma1)) + (mu1-mu2) %*% solve(Sigma2) %*% (mu1-mu2) ) / 2)
}
