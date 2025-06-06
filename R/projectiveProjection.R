#' Projective projection
#'
#' @param Y High dimensional data to project
#' @param p Dimensionality to project to
#' @param spherize Spherize dataset or not
#' 
#' @import corpcor
#'
#' @return matrix with coordinates of projected points.
projectiveProjection <- function(Y, p, spherize=F) {
  L <- nrow(Y)
  N <- ncol(Y)
  
  Ymean <- apply(Y, 1, mean)
  ym <- matrix(Ymean, ncol=1)
  Y <- Y - Ymean
  svdObj <- fast.svd(Y)
  Up <- svdObj$u[, 1:(p-1)]
  proj <- Up %*% t(Up)
  D <- svdObj$d[1:(p-1)]
  
  Y <- proj %*% Y
  Y <- Y + Ymean
  YmeanOrtho <- ym - proj %*% ym
  Up <- cbind(Up, YmeanOrtho / (sqrt(sum(YmeanOrtho ^ 2))))
  singValues <- D
  lamSphe <- 1e-8
  
  Y <- t(Up) %*% Y
  
  ## spherizing
  if (spherize) {
    Y <- Up %*% Y
    Y <- Y - Ymean
    C <- diag(1 / sqrt(D + lamSphe))
    IC <- solve(C)
    Y <- C %*% t(Up[, 1:(p-1)]) %*% Y
    Y <- rbind(Y, 1)
    Y <- Y / sqrt(p)
  }
  return(Y)
}


#' Get projective projection
#'
#' @param Y High dimensional data to project
#' @param p Dimensionality to project to
#' @param spherize Spherize dataset or not
#' 
#' @import corpcor
#'
#' @return projection
getProjectiveProjection <- function(Y, p, spherize=F) {
  L <- nrow(Y)
  N <- ncol(Y)
  
  Ymean <- apply(Y, 1, mean)
  ym <- matrix(Ymean, ncol=1)
  Y <- Y - Ymean
  svdObj <- fast.svd(Y)
  Up <- svdObj$u[, 1:(p-1)]
  proj <- Up %*% t(Up)
  D <- svdObj$d[1:(p-1)]
  
  Y <- proj %*% Y
  Y <- Y + Ymean
  YmeanOrtho <- ym - proj %*% ym
  Up <- cbind(Up, YmeanOrtho / (sqrt(sum(YmeanOrtho ^ 2))))

  return(Up)
}