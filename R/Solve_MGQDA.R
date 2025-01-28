#' Solve_MGQDA Function
#'
#' Apply the block-coordinate descent algorithm and Newton's method to obtain the projection matrix V.
#'
#' @param X An n by p matrix representing the samples, where n is the number of samples, and p is the number of features.
#'
#' @param Y A numeric vector of length n, representing the group to which each sample belongs.
#'
#' @param alpha The penalty coefficient alpha.
#'
#' @param lambda The penalty coefficient lambda.
#'
#' @param N The number of iterations.
#'
#' @param eps The tolerance level for iterative convergence.
#'
#' @param maxiter The maximum number of iterations allowed.
#'
#' @return Projection matrix V, features number nfeature, iteration times iter.
#'
#' @import MASS
#'
#' @export


Solve_MGQDA<- function(X, Y, alpha, lambda, V = NULL, eps = 1e-4, N = 0, maxiter = 10000){
  if (any(is.na(X))|any(is.na(Y))){
    stop("Missing values are not allowed")
  }
  p <- ncol(X)
  n <- length(Y)
  G <- max(Y)
  if (nrow(X) != n){
    stop("Dimensions of X and Y don't match!")
  }

  if (any(apply(X,2,sd) < 1e-13)){
    stop(paste("Some features have standard deviation less than 1e-13!", sep = ""))
  }

  if (is.null(V)){
    V = matrix(0, p, G*(G-1))
  }

  if(G > 2){
    YT = list()
    for(i in 1:G)
      YT[[i]] = which(Y == i)
    Yt <- rep(0,G)
    for(i in 1:G)
      Yt[i] = length(YT[[i]])
    XT = list()
    for(i in 1:G)
      XT[[i]] = X[YT[[i]],]
    XTMeans = list()
    for(i in 1:G)
      XTMeans[[i]] = colMeans(XT[[i]])
    XMeans = colMeans(X)
    B = matrix(rep(0, p*p), p, p)
    for(i in 1:G)
      B = B + (Yt[i]/n)*tcrossprod(XTMeans[[i]]-XMeans)
    D = matrix(rep(0,p*(G-1)), nrow = p, ncol = (G-1))
    if(G == 2){
      D = (sqrt(Yt[2])*Yt[1]*(XTMeans[[1]]-XTMeans[[2]]))/(n*sqrt(Yt[1]))
    }
    if(G >= 3){
      deltaXmean <- matrix(rep(0,p*(G-1)), nrow = p, ncol = (G-1))
      for(i in 1:(G-1))
        deltaXmean[,i] = Yt[i]*(XTMeans[[i]])
      for(i in 1:(G-1)){
        if(i == 1){
          D[,i] = (sqrt(Yt[2])*Yt[1]*(XTMeans[[1]]-XTMeans[[2]]))/(sqrt(n)*sqrt(Yt[1])*sqrt(Yt[1]+Yt[2]))
        }
        if(i >= 2){
          D[,i] = (sqrt(Yt[i+1])*((rowSums(deltaXmean[,1:i])) - sum(Yt[1:i])*XTMeans[[i+1]]))/(sqrt(n)*sqrt(sum(Yt[1:i])*sum(Yt[1:(i+1)])))
        }
      }
    }
    Sig = crossprod(X-matrix(rep(colMeans(X),n),n,p,byrow = TRUE))
    S = list()
    for(i in 1:G)
      S[[i]] = crossprod(XT[[i]]-matrix(rep(colMeans(XT[[i]]),Yt[i]),Yt[i],p,byrow = TRUE))/Yt[i]
    S_max = matrix(0, p, p*G)
    for(i in 1:G)
      S_max[,((i-1)*p+1):(i*p)] = S[[i]]
    out = updateMatrix(as.matrix(S_max), as.matrix(B), as.matrix(D), as.matrix(V),as.double(lambda), as.double(alpha), as.double(eps), maxiter = maxiter, p = p, G = G)
    V = out$V
    nfeature = sum(rowSums(abs(V))>0)
    iter = out$iter
    return(list(V = V, nfeature = nfeature, iter = iter))
  }
  if(G == 2){

    YT = list()
    YT[[1]] = which(Y == 1)
    YT[[2]] = which(Y == 2)

    Yt <- rep(0,2)
    Yt[1] = length(YT[[1]])
    Yt[2] = length(YT[[2]])

    XT = list()
    XT[[1]] = X[YT[[1]],]
    XT[[2]] = X[YT[[2]],]

    XTMeans = list()
    XTMeans[[1]] = colMeans(XT[[1]])
    XTMeans[[2]] = colMeans(XT[[2]])

    XMeans = colMeans(X)

    B = matrix(rep(0, p*p), p, p)
    B = B + (Yt[1]/n)*tcrossprod(XTMeans[[1]]-XMeans) + (Yt[2]/n)*tcrossprod(XTMeans[[2]]-XMeans)

    D = rep(0,p*(G-1))
    D = (sqrt(Yt[2])*Yt[1]*(XTMeans[[1]]-XTMeans[[2]]))/(n*sqrt(Yt[1]))

    Sig = crossprod(X-matrix(rep(colMeans(X),n),n,p,byrow = TRUE))

    S = list()
    for(i in 1:G)
      S[[i]] = crossprod(XT[[i]]-matrix(rep(colMeans(XT[[i]]),Yt[i]),Yt[i],p,byrow = TRUE))
    S_max = matrix(0, p, p*G)
    for(i in 1:G)
      S_max[,((i-1)*p+1):(i*p)] = S[[i]]
    out = updateMatrix(as.matrix(S_max), as.matrix(B), as.vector(D), as.matrix(V),as.double(lambda), as.double(alpha), as.double(eps), maxiter = maxiter, p = p, G = G)
    V = out$V
    nfeature = sum(rowSums(abs(V))>0)
    iter = out$iter
    return(list(V = V, nfeature = nfeature, iter = iter))
  }
}
