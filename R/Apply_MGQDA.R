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
