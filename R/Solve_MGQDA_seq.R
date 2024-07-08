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
