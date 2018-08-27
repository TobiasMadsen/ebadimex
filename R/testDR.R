#'
#' @export
testDR <- function(mat, grouping, prior){
  ###
  # Compute projections and SSD matrices
  ###
  T1 <- cbind(1, grouping)
  P1 <- T1 %*% solve(t(T1) %*% T1) %*% t(T1)
  SSD1 <- (t(mat) %*% mat - t(mat) %*% P1 %*% mat)
  
  T2 <- rep(1, length(grouping))
  P2 <- T2 %*% solve(t(T2) %*% T2) %*% t(T2)
  SSD2 <- t(mat) %*% (P1 - P2) %*% mat
  
  ###
  # Perform tests
  ###
  df1 <- nrow(mat) - 2
  p <- ncol(mat)
  
  Q2n <- det(SSD1 + prior$nu0 * prior$Lambda) / det(SSD1 + SSD2 + prior$nu0 * prior$Lambda)
  f1 <- nrow(mat) - 2 + prior$nu0 - p + 1
  f2 <- p 
  fStat <- (1 - Q2n) / Q2n *  f1 / f2
    
  pObs <- pf(fStat, df1 = f2, df2 = f1, lower.tail = F, log.p = T)
  
  pObs
}