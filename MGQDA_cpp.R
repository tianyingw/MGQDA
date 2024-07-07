### MGQDA by Rcpp
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("update_matrix_MGQDA.cpp")

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
Solve_MGQDA_seq <- function(X, Y, alpha, lambda_seq, eps = 1e-4, maxiter = 10000, feature_max = nrow(X)){
  p = ncol(X)
  G = max(Y)
  # Number of lambda in given data (some of them may needless to compute because too small cannot guarantee the sparsity, function will terminate prematurely)
  n_lambda = length(lambda_seq)
  # This function V is a bind function: fromn left to right consist n_lambda different matrix V, corresponding to each different lambda
  V = matrix(NA, nrow = p, ncol = (n_lambda*G*(G-1)))
  # Starting value for first iteration
  V0 = matrix(0, nrow = p, ncol = (G*(G-1)))
  nfeature_vec = rep(NA, n_lambda)
  for(i in 1:n_lambda){
    out = Solve_MGQDA(X, Y, alpha = alpha, lambda = lambda_seq[i], V = V0)
    # Previous update V0 as a warm start
    V0 = out$V
    # Update a new p x G(G-1) matrix in whole matrix
    V[,c(((i-1)*G*(G-1)+1):(i*G*(G-1)))] = out$V
    nfeature_vec[i] = out$nfeature
    # Terminate the function if sparsity cannot hold
    if(out$nfeature >= feature_max){
      return(list(V = V[,c(1:(i*G*(G-1)))], lambda_seq = lambda_seq[1:i], nfeature_vec = nfeature_vec[1:i]))
    }
  }
  return(list(V = V, lambda_seq = lambda_seq, nfeature_vec = nfeature_vec))
}
Pred_QDA <- function(Xtrain, Ytrain, Xtest, Ytest = NULL, V = NULL, prior = FALSE){
  if(!is.null(V)){
    p <- ncol(Xtrain)
    ntrain <- length(Ytrain)
    G <- max(Ytrain)
    ntest = nrow(Xtest)
    ## Use svd to obtain singular value of matrix V : V = Usvd %*% diag(d) %*%t(Vsvd)
    Vsvd = base::svd(V)$d
    Vsvdu = base::svd(V)$u
    Vsvdv = base::svd(V)$v
    if(min(Vsvd) < 1e-6){
      Vsvd = Vsvd[-which(Vsvd < 1e-6)]
    }
    m = length(Vsvd)
    ## Three different cases : m=1; m=2~G(G-1)-1; m=G(G-1).
    ## For general case: m=2~G(G-1), can directly apply the package.
    if(m >= 2){
      V = Vsvdu[,c(1:m)] %*% diag(Vsvd) %*% t(Vsvdv[c(1:m), c(1:m)])
      ## Eliminate elements small enough
      for(i in 1:p){
        for(j in 1:m){
          if(abs(V[i,j]) < 1e-13){
            V[i,j] = 0
          }
        } 
      }
      Xtrain = Xtrain %*% V
      Xtest = Xtest %*% V
      if(is.null(dim(Xtrain))) stop("'Xtrain' is not a matrix")
      if(is.null(dim(Xtest))) stop("'Xtest' is not a matrix")
      Xtrain <- as.matrix(Xtrain)
      Xtest <- as.matrix(Xtest)
      if(any(!is.finite(Xtrain)))
        stop("infinite, NA or NaN values in 'Xtrain'")
      if(any(!is.finite(Xtest)))
        stop("infinite, NA or NaN values in 'Xtest'")
      if(ncol(Xtrain) != ncol(Xtest))
        stop("incorrect dimension between Xtrain and Xtest")
      if(nrow(Xtrain) != length(Ytrain))
        stop("counts of training set doesn't match")
      p <- ncol(Xtrain)
      ntrain <- nrow(Xtrain)
      G <- max(Ytrain)
      Yt <- rep(0,p)
      for(i in 1:p){
        Yt[i] = sum(Ytrain == i)
      }
      Xt = list()
      for(i in 1:G){
        Xt[[i]] = Xtrain[which(Ytrain == i),]
      }
      Sg = list()
      for(i in 1:G){
        Sg[[i]] = crossprod(t((t(Xt[[i]])-colMeans(Xt[[i]]))))/(sum(Ytrain == i)-1)
      }
      Xtmeans = list()
      for(i in 1:G){
        Xtmeans[[i]] = colMeans(Xt[[i]])
      }
      Sginv = list()
      for(i in 1:G){
        Sginv[[i]] = ginv(Sg[[i]])
      }
      LnSg = list()
      for(i in 1:G){
        LnSg[[i]] = log(abs(det(Sg[[i]]))+1e-15)
      }
      ntest = nrow(Xtest)
      testmat <- matrix(1, ntest, G)
      if(prior == FALSE){
        for(i in 1:ntest){
          for(j in 1:G){
            testmat[i,j] = (t(Xtest[i,]) - Xtmeans[[j]])%*%Sginv[[j]]%*%t(t(Xtest[i,]) - Xtmeans[[j]]) + LnSg[[j]] - 2*sum(Ytrain == j)/ntrain
          }
        }
      }else{
        for(i in 1:ntest){
          for(j in 1:G){
            testmat[i,j] = (t(Xtest[i,]) - Xtmeans[[j]])%*%Sginv[[j]]%*%t(t(Xtest[i,]) - Xtmeans[[j]]) + LnSg[[j]]
          }
        }
      }
      pred_vec = rep(0,ntest)
      for(i in 1:ntest){
        pred_vec[i] = which.min(testmat[i,])
      }
      if(!is.null(Ytest)){
        if(length(Ytest) == ntest & max(Ytest) == G){
          errorcounts = sum(pred_vec != Ytest)
          errorrate = errorcounts/ntest
          return(list(pred_vec = pred_vec, errorcounts = errorcounts, errorrate = errorrate))
        }else{
          return(pred_vec = pred_vec)
        }
      }else{
        return(pred_vec = pred_vec)
      }
      
    }
    if (m == 1){
      V = Vsvdu[,1]*Vsvd*Vsvdv[1,1]
      for(i in 1:p){
        if(abs(V[i]) < 1e-13){
          V[i] = 0
        }
      }
      XtrainProj = Xtrain %*% V
      XtestProj = Xtest %*% V
      Omega = list()
      for(i in 1:G)
        Omega[[i]] = 1 / stats::var(XtrainProj[Ytrain == i, ])
      mu = list()
      for(i in 1:G)
        mu[[i]] = mean(XtrainProj[Ytrain == i, ])
      ntest = nrow(Xtest)
      Predict_mat = matrix(1, ntest, G)
      for(i in 1:G){
        Predict_mat[,i] = (XtestProj-mu[[i]])^2*Omega[[i]] - log(Omega[[i]])
      }
      if(prior == TRUE){
        for(i in 1:G){
          Predict_mat[,i] = Predict_mat[,i] - 2*log(length(Ytrain[which(Ytrain == i)])/ntrain)
        }
      }
      pred_vec = rep(1,ntest)
      for(i in 1:ntest){
        pred_vec[i] = which.min(Predict_mat[i,])
      }
      if(!is.null(Ytest)){
        if(length(Ytest) == ntest & max(Ytest) == G){
          errorcounts = sum(pred_vec != Ytest)
          errorrate = errorcounts/ntest
          return(list(pred_vec = pred_vec, errorcounts = errorcounts, errorrate = errorrate))
        }else{
          return(pred_vec = pred_vec)
        }
      }else{
        return(pred_vec = pred_vec)
      }
    }
    if (m == 0){
      pred_vec = rep(1,ntest)
      pred_vec = floor(runif(nrow(Xtest), 1, (G+1)))
      if(!is.null(Ytest)){
        if(length(Ytest) == ntest & max(Ytest) == G){
          errorcounts = sum(pred_vec != Ytest)
          errorrate = errorcounts/ntest
          return(list(pred_vec = pred_vec, errorcounts = errorcounts, errorrate = errorrate))
        }else{
          return(pred_vec = pred_vec)
        }
      }else{
        return(pred_vec = pred_vec)
      }
    }
  }else{
    if(is.null(dim(Xtrain))) stop("'Xtrain' is not a matrix")
    if(is.null(dim(Xtest))) stop("'Xtest' is not a matrix")
    Xtrain <- as.matrix(Xtrain)
    Xtest <- as.matrix(Xtest)
    if(any(!is.finite(Xtrain)))
      stop("infinite, NA or NaN values in 'Xtrain'")
    if(any(!is.finite(Xtest)))
      stop("infinite, NA or NaN values in 'Xtest'")
    if(ncol(Xtrain) != ncol(Xtest))
      stop("incorrect dimension between Xtrain and Xtest")
    if(nrow(Xtrain) != length(Ytrain))
      stop("counts of training set doesn't match")
    p <- ncol(Xtrain)
    ntrain <- nrow(Xtrain)
    G <- max(Ytrain)
    Yt <- rep(0,p)
    for(i in 1:p){
      Yt[i] = sum(Ytrain == i)
    }
    Xt = list()
    for(i in 1:G){
      Xt[[i]] = Xtrain[which(Ytrain == i),]
    }
    Sg = list()
    for(i in 1:G){
      Sg[[i]] = crossprod(t((t(Xt[[i]])-colMeans(Xt[[i]]))))/(sum(Ytrain == i)-1)
    }
    Xtmeans = list()
    for(i in 1:G){
      Xtmeans[[i]] = colMeans(Xt[[i]])
    }
    Sginv = list()
    for(i in 1:G){
      Sginv[[i]] = ginv(Sg[[i]])
    }
    LnSg = list()
    for(i in 1:G){
      LnSg[[i]] = log(abs(det(Sg[[i]]))+1e-15)
    }
    ntest = nrow(Xtest)
    testmat <- matrix(1, ntest, G)
    if(prior == FALSE){
      for(i in 1:ntest){
        for(j in 1:G){
          testmat[i,j] = (t(Xtest[i,]) - Xtmeans[[j]])%*%Sginv[[j]]%*%t(t(Xtest[i,]) - Xtmeans[[j]]) + LnSg[[j]] - 2*sum(Ytrain == j)/ntrain
        }
      }
    }else{
      for(i in 1:ntest){
        for(j in 1:G){
          testmat[i,j] = (t(Xtest[i,]) - Xtmeans[[j]])%*%Sginv[[j]]%*%t(t(Xtest[i,]) - Xtmeans[[j]]) + LnSg[[j]]
        }
      }
    }
    pred_vec = rep(0,ntest)
    for(i in 1:ntest){
      pred_vec[i] = which.min(testmat[i,])
    }
    if(!is.null(Ytest)){
      if(length(Ytest) == ntest & max(Ytest) == G){
        errorcounts = sum(pred_vec != Ytest)
        errorrate = errorcounts/ntest
        return(list(pred_vec = pred_vec, errorcounts = errorcounts, errorrate = errorrate))
      }else{
        return(pred_vec = pred_vec)
      }
    }else{
      return(pred_vec = pred_vec)
    }
  }
}
CV_MGQDA <- function(X, Y, alpha, lambda_seq, nfolds = 5, eps = 1e-4, maxiter = 1000, myseed = 1001, prior = TRUE){
  
  n = length(Y)
  n_lambda = length(lambda_seq)
  G = max(Y)
  
  error_mat = matrix(1, nfolds, n_lambda)
  error_sum = matrix(NA, nfolds, n_lambda)
  
  ####random set split the whole data set into k folds corresponding to ng
  set.seed(6700417)
  id = rep(NA, n)
  for(i in 1:G){
    id[Y == i]= sample(rep(seq_len(nfolds), length.out = sum(Y == i)))
  }
  
  
  nfeature_mat = matrix(NA, nfolds, n_lambda)
  for (nf in 1: nfolds){
    Xtrain = X[id != nf, ]
    Ytrain = Y[id != nf]
    Xtest = X[id == nf, ]
    Ytest = Y[id == nf]
    ## For this fold: different lambda leads to different V and nfeature_vec
    fit_tmp = Solve_MGQDA_seq(X = Xtrain, Y = Ytrain, alpha = alpha, lambda_seq = lambda_seq, eps = 1e-4, maxiter = 10000, feature_max = n)
    nfeature_mat[nf, 1:length(fit_tmp$nfeature_vec)] = fit_tmp$nfeature_vec
    l = length(fit_tmp$lambda_seq)
    ## Obtain the V we need in cross-validation
    V = fit_tmp$V
    Vc = list()
    if(l == 1){
      Vc[[1]] = V
    }else{
      for(i in 1:l){
        Vc[[i]] = V[,c(((i-1)*G*(G-1)+1):(i*G*(G-1)))]
      }
    }
    for(j in 1:length(fit_tmp$lambda_seq)){
      Ypred <- Pred_QDA(Xtrain = Xtrain, Ytrain = Ytrain, Xtest = Xtest, V = Vc[[j]])
      error_mat[nf, j] = sum(Ypred != Ytest) / length(Ytest)
      error_sum[nf, j] = sum(Ypred != Ytest) 
    }
  }
  cvm = colSums(error_sum)/length(Y)
  index = which.min(cvm)
  lambda_min = lambda_seq[index]
  sd_vec = rep(0, n_lambda)
  for(j in 1: n_lambda){
    sd_vec[j] = sd(error_mat[,j])
  }
  cv_sd = cvm[index] + sd_vec[index]
  sd_index = which(cvm <= cv_sd)[1]
  lambda_sd = lambda_seq[sd_index]
  
  
  return (list(lambda_min = lambda_min, lambda_sd = lambda_sd, error_mat = error_mat, cv_sd = cv_sd, sd_vec = sd_vec))
}
Apply_MGQDA <- function(Xtrain, Ytrain, Xtest, Ytest = NULL, alpha = NULL, lambda_seq = NULL, n_lambda = 25,  maxmin_ratio = 0.2, nfolds = 5, eps = 1e-5, maxiter = 10000, myseed = 1001, prior = TRUE){
  if (!is.null(alpha)) {
    nl1 = length(alpha)
    if (nl1 > 1){
      warning(paste("The parameters you provided are not compliant. Alpha will be set as the first element."))
      alpha = alpha[1]
    }else{
      alpha = alpha
    }
  }else {
    alpha = 0.5
  }
  ## Find proper lambda_max by bisection method
  l_max = 1
  G = max(Ytrain)
  if(Solve_MGQDA(X = Xtrain, Y = Ytrain, alpha = alpha, lambda = l_max)$nfeature > (G*(G-1))){
    repeat{
      l_max = 1.25*l_max
      if(Solve_MGQDA(X = Xtrain, Y = Ytrain, alpha = alpha, lambda = l_max)$nfeature <= (G*(G-1)))break
    }
    l_max = 0.8*l_max
  }
  if(Solve_MGQDA(X = Xtrain, Y = Ytrain, alpha = alpha, lambda = l_max)$nfeature <= (G*(G-1))){
    repeat{
      l_max = 0.8*l_max
      if(Solve_MGQDA(X = Xtrain, Y = Ytrain, alpha = alpha, lambda = l_max)$nfeature > (G*(G-1)))break
    }
  }
  ## Generate lambda sequence
  if (!is.null(lambda_seq)) {
    nl = length(lambda_seq)
    if (nl < 1){
      warning(paste("There is no qualified lambda value. New values will be generated automatically. n_lambda will be set as.", n_lambda,sep = " "))
      lambda_seq = exp(seq(log(l_max * maxmin_ratio), log(l_max), length.out = n_lambda))
    }else{
      n_lambda = nl
    }
  }else {
    lambda_seq = exp(seq(log(l_max * maxmin_ratio), log(l_max), length.out = n_lambda))
  }
  lambda_seq = sort(lambda_seq, decreasing = FALSE)
  ## Use CV to select the tuning parameter
  out.cv = CV_MGQDA(X = Xtrain, Y = Ytrain, alpha = alpha, lambda_seq = lambda_seq, nfolds = nfolds, eps = eps, maxiter = maxiter, myseed = myseed, prior = prior)
  out.proj = Solve_MGQDA(X = Xtrain, Y = Ytrain, alpha = alpha, lambda = out.cv$lambda_min, eps = eps, maxiter = maxiter)
  V = out.proj$V
  if (out.proj$nfeature > 0){
    Ypred = Pred_QDA(Xtrain = Xtrain, Ytrain = Ytrain, Xtest = Xtest, V = V, prior = prior)
    if (!is.null(Ytest))
    {
      error = sum(Ypred != Ytest) / length(Ytest)
      errorG = rep(1, max(Ytest))
      for(g in 1:(max(Ytest))){
        errorG[g] = sum(Ytest[which(Ytest==g)]!=Ypred[which(Ytest==g)])/length(Ytest[which(Ytest==g)])
      }
      return(list(error = error, errorGroup = errorG, features = out.proj$nfeature, features_id = c(1:ncol(Xtrain))[rowSums(abs(V)) != 0],  lambda = out.cv$lambda_min, V = V))
    }else{
      return(list(Ypred = Ypred, features = out.proj$nfeature, features_id = c(1:ncol(Xtrain))[rowSums(abs(V)) != 0],  lambda = out.cv$lambda_min, V = V))
    }
  }else{
    return(list(error = 1, features = 0, features_id = NA))
  }
}















s = Sys.time(); MGQDARES1 <- MGQDAsimulation1(p = 100,ng = 100,rho = 0.8, b = 30, rp = 3); t = Sys.time(); t-s