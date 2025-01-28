#' CV_MGQDA Function
#'
#' Select the optimal penalty coefficient through cross-validation.
#'
#' @param X An n by p matrix representing the samples, where n is the number of samples, and p is the number of features.
#'
#' @param Y A numeric vector of length n, representing the group to which each sample belongs.
#'
#' @param alpha The penalty coefficient alpha.
#'
#' @param lambda_seq A sequence of lambda values for cross-validation, used to select an appropriate lambda.
#'
#' @param nfolds The number of folds to use for cross-validation.
#'
#' @param eps The tolerance level for iterative convergence.
#'
#' @param maxiter The maximum number of iterations allowed.
#'
#' @param myseed The random seed used for reproducibility.
#'
#' @param prior Logical; if TRUE, prior probabilities are used.
#'
#' CV_MGQDA(X, Y, alpha, lambda_seq, nfolds = 5, eps = 1e-4, maxiter = 1000, myseed = 1001, prior = TRUE)
#' @export


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

