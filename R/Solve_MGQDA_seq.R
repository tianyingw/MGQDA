#' Solve_MGQDA_seq Function
#'
#' Apply a sequence of lambda values to the Solve_MGQDA function to select the appropriate penalty coefficient lambda.
#'
#' @param X An n by p matrix representing the samples, where n is the number of samples, and p is the number of features.
#'
#' @param Y A numeric vector of length n, representing the group to which each sample belongs.
#'
#' @param alpha The penalty coefficient alpha.
#'
#' @param lambda_seq A sequence of lambda values, used to select an appropriate lambda through cross-validation.
#'
#' @param eps The tolerance level for iterative convergence.
#'
#' @param feature_max The maximum number of features to be selected.
#'
#' @param maxiter The maximum number of iterations allowed.
#'
#' Solve_MGQDA_seq(X, Y, alpha, lambda_seq, eps = 1e-4, maxiter = 10000, feature_max = nrow(X))
#' @export

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
