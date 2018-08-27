#' @param varEst matrix of dimension N x k containing variance estimates for each measurements type and each unit. 
#' N is the number of units and k is the number of measurements for each unit.
#' @param d degrees of freedom in conditional Wishart distribution for variance/covariance estimates
fitInverseChiSquare <- function(varEst, d){
  # Cast to matrix
  varEst <- as.matrix(varEst)
  
  # Using notation similar to (Smyth 2004)
  z <- log(varEst)  
  e <- z - digamma(d/2) + log(d/2)
  
  # Estimate Degrees of freedom
  solve_trigamma_nu0_half <- mean( apply(e, MARGIN = 2, FUN = var)  ) - trigamma(d/2)
  f <- (ncol(varEst)-1)
  nu0 <- uniroot(function(x){trigamma(x) - solve_trigamma_nu0_half}, interval = c(1e-9, 1000))$root*2
  nu0 <- nu0 + f
  
  # Estimate locations
  s0 <- exp( (colMeans(e) + digamma( (nu0-f)/2) - log( (nu0-f) / 2)) ) * (nu0-f) / nu0
  
  return(list(s0 = s0, nu0 = nu0))
}

#' Estimate Prior for Covariance Matrices
#'
#'
#' @export
fitInverseWishart <- function(covMats, d){
  # Estimate degrees of freedom and diagonal entries of Lambda_0
  varEst <- t(sapply(covMats, FUN = diag))
  fit <- fitInverseChiSquare(varEst, d)
 
  # Estimate off-diagonal entries using moment matching to 
  # fisher Z transformation (atanh) to correlation estimates
  
  #corMats <- lapply(covMats, FUN = cov2cor)
  #zCorMats <- lapply(corMats, FUN = atanh)
  #meanZCorMats <- Reduce('+', zCorMats) / length(zCorMats)
  #corMatEst <- tanh(meanZCorMats)
  corMatEst <- cov2cor( Reduce('+', covMats) / length(covMats) * (fit$nu0 - 2 - 1) / fit$nu0 ) 

  return(list(Lambda =  corMatEst * sqrt(fit$s0) * rep(sqrt(fit$s0), each = nrow(corMatEst)),
              nu0 = fit$nu0 ) )
}

#' Find prior for expression data
#' 
#' @param exprMethList List of nx3 matrices with expression measurements in the first column and
#' methylation measurements in second and third column
#'
#' @import ggplot2
#' 
#' @export
exprPrior <- function(exprMethList, bins = 10, plot = TRUE){
  # Fit some class of distribution parameterised by mean and variance to each gene
  m <- sapply(exprMethList, function(x){mean(x[,1])})
  v <- sapply(exprMethList, function(x){var(x[,1])})
  numberOfSamples <- nrow(exprMethList[[1]])
  d <- numberOfSamples-1
  
  # Fit inverse gamma for different bins
  quan_m <- seq(min(m), max(m), length.out = bins + 1)
  quan_m <- quantile(m, probs = seq(0,1,length.out = bins +1))
  quan_s_fit <- list()
  param_s_fit <- list()
  for(i in 1:bins){
    idx <- quan_m[i] < m & m < quan_m[i+1] 
    fit <- fitInverseChiSquare(v[idx], d)
    
    param_s_fit[[i]] <- c(mean = unname((quan_m[i]+quan_m[i+1])/2), s0 = fit[['s0']], d0 = fit[['nu0']])
    quan_s_fit[[i]] <- log(fit[['s0']]) + log(qf(p = seq(0.05,0.95,0.05), df1 = d, df2 = fit[['nu0']]))
  }
  
  if (plot) {
    # Plot 
    plot_df <- data.frame(cbind(do.call(rbind, quan_s_fit))[c(1,1:bins,bins),], m = c(min(m), (quan_m[2:(bins+1)] + quan_m[1:bins])/2, max(m)) )
    names(plot_df)[1:19] <- paste0('q', seq(5,95,5))
    p <- ggplot(data.frame(x = m, y = log(v)), aes(x = x, y = y)) + geom_point(alpha = 0.1) +
      geom_ribbon(data = plot_df, aes(x = m, ymin = q25, ymax = q75, y = NULL), alpha = 0.4, fill = "aquamarine4") +
      geom_line(data = plot_df, aes(x = m, y = q50), colour = "dodgerblue", linetype = "dashed") +
      geom_line(data = plot_df, aes(x = m, y = q25), alpha = 0.3) +
      geom_line(data = plot_df, aes(x = m, y = q75), alpha = 0.3) +
      ylab("Log variance") +
      xlab("Mean log expression") +
      theme_minimal()
    plot(p)
  }
  
  # Return
  ret <- do.call(rbind,param_s_fit)
  colnames(ret) <- c("expr", "s0", "d0")
  return(ret)
}

#' Find prior for methylation data
#' @param exprMethList List of nx3 matrices with expression measurements in the first column and
#' methylation measurements in second and third column
#'
#' @export
methPrior <- function(exprMethList, grouping){
  # Estimate SSD Matrices
  d <- length(grouping) - 3
  
  covMatList <- lapply(1:length(exprMethList), function(i){
    meth <- exprMethList[[i]][,2:3]
    expr <- exprMethList[[i]][,1]
    fitMeth(meth, expr, grouping )$SSD1 / d
  })
  
 return(fitInverseWishart(covMatList,d))
}