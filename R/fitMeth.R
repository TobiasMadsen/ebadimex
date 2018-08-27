#' @import dplyr
#' @importFrom reshape2 melt
#' @export
methToLongFormat <- function(meth, expr, grouping = NULL){
  if(is.null(grouping)){
    melt(meth, varnames = c("sample", "probe")) %>% mutate(type = ifelse(str_detect(as.character(probe), "pr"), "pr", "gb")) %>%
      mutate(expr = expr[ (row_number()-1) %% nrow(meth) + 1]) %>%
      mutate(pr = as.numeric(type == "pr"), gb = as.numeric(type == "gb")) %>%
      mutate(type = as.factor(type)) %>% return
  } else {
    melt(meth, varnames = c("sample", "probe")) %>% mutate(type = ifelse(str_detect(as.character(probe), "pr"), "pr", "gb")) %>%
      mutate(expr = expr[ (row_number()-1) %% nrow(meth) + 1]) %>% 
      mutate(group = grouping[(row_number()-1) %% nrow(meth) + 1]) %>%
      mutate(pr = as.numeric(type == "pr"), gb = as.numeric(type == "gb")) %>%
      mutate(diff_pr = (type == "pr") * group) %>%
      mutate(diff_gb = (type == "gb") * group) %>%
      mutate(type = as.factor(type)) %>% return
  }
}


#' Fit expression data Linear Normal Model
#' @param meth matrix with met
#' @param meth matrix with methylation values
#' @param expr vector with expression values
#' @param grouping logical vector indicating group relationship
#' 
#' 
#' @importFrom stringr str_detect
#' 
#' @export
fitMeth <- function(meth, expr, grouping, win = F, prior = NULL){
  if(is.null(prior)){
    prior <- list(Lambda = matrix(0,2,2),
                  nu0 = 0)
  }
  
  # Determine probe status: Promoter or Genebody
  N <- nrow(meth)
  n1 <- sum(grouping)
  n2 <- sum(!grouping)
  
  
  # Winsorize
  if(win){
    meth[,1] <- winsorize(meth[,1], g = grouping, k_win = 2.5, boost = TRUE)$x_star
    meth[,2] <- winsorize(meth[,2], g = grouping, k_win = 2.5, boost = TRUE)$x_star
  }
  
  # Design matrices
  Tm1 <- cbind(expr, 1, ifelse(grouping, 1, 0))    
  Tm2 <- cbind(expr, 1)
  
  # Learn parameters
  alpha1 <- solve(t(Tm1) %*% Tm1 ) %*% t(Tm1) %*% meth
  alpha2 <- solve(t(Tm2) %*% Tm2 ) %*% t(Tm2) %*% meth
  
  # Covariance matrix
  X1 <- meth[grouping,]
  X2 <- meth[!grouping,]
  SSDg1 <- t(X1) %*% X1 - t(X1) %*% Tm1[grouping,]  %*% alpha1
  SSDg2 <- t(X2) %*% X2 - t(X2) %*% Tm1[!grouping,]  %*% alpha1
  SSD1 <- SSDg1 + SSDg2

  SSD02 <- t(meth) %*% meth - t(meth) %*% Tm2 %*% alpha2
  
  # Regularize estimates
  SSD1_reg  <- SSD1 + prior$nu0 * prior$Lambda
  SSDg1_reg <- SSDg1 + prior$nu0 * prior$Lambda
  SSDg2_reg <- SSDg2 + prior$nu0 * prior$Lambda
  
  # Ensure symmetry
  SSD1[1,2] <- SSD1[2,1]
  SSDg1[1,2] <- SSDg1[2,1]
  SSDg2[1,2] <- SSDg2[2,1]
  SSD1_reg[1,2] <- SSD1_reg[2,1]
  SSDg1_reg[1,2] <- SSDg1_reg[2,1]
  SSDg2_reg[1,2] <- SSDg2_reg[2,1]
  
  # Return parameters and SSD Matrix
  param <- list()
  param$a_pr <- alpha1[1,1]
  param$a_gb <- alpha1[1,2]
  param$alpha <- alpha1
  param$Sigma <- SSD1 / (N-3)
  param$Sigma_g1 <- SSDg1 / (n1-2)
  param$Sigma_g2 <- SSDg2 / (n2-2)
  param$Sigma_reg   <- SSD1_reg / (N+prior$nu0-3)
  param$Sigma_g1_reg  <- SSDg1_reg / (n1+prior$nu0-2)
  param$Sigma_g2_reg  <- SSDg2_reg / (n2+prior$nu0-2)
  
  param$diff_x <- alpha1[3,1]
  param$diff_y <- alpha1[3,2]
  
  # Loglik 
  loglik <- -N*log(2*pi)-N/2*det(SSD1/N)-N
  
  return(list(param = param,
              SSD1 = SSD1,
              SSD2 = SSD02 - SSD1 ,
              loglik = loglik,
              SSDg1 = SSDg1,
              SSDg2 = SSDg2,
              Lambda_0 = prior$Lambda,
              nu0 = prior$nu0))
}
