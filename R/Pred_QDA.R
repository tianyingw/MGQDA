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
