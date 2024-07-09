# Description of package MGQDA

The main funcitons are CV_MGQDA (cross-validation), Pred_QDA (prediction) and Solve_MGQDA (apply the rule of selected tuning parameter).

### Simple Example

```R
library(MGQDA)

### Example 1
# generate training data
n <- 100
p <- 200
G <- 3
Ytrain <- rep(1:G, each = n)
set.seed(123)
Xtrain <- matrix(rnorm(p * n * G), n * G, p)

# find projection matrix V with selected tuning parameter alpha = 0.1 and lambda = 0.5:
V <- Solve_MGQDA(Xtrain, Ytrain, alpha = 0.1, lambda = 0.5)
nfeature = sum(rowSums(V) != 0)

# generate test data
m <- 20
set.seed(321)
Xtest <- matrix(rnorm(p * m), m, p)

# perform classification
Ytest <- Pred_QDA(Xtrain, Ytrain, Xtest, V)

# select tuning parameter automatically
out = Apply_MGQDA(Xtrain, Ytrain, Xtest, Ytest = NULL, alpha = NULL, lambda_seq = NULL, n_lambda = 25,  maxmin_ratio = 0.2, nfolds = 5, eps = 1e-5, maxiter = 10000, myseed = 1001, prior = TRUE)
# total error rate
error = out$error
# errorrate for each group
errorGroup = out$errorGroup
# number of selected features and features id
features = out$features
features_id = out$features_id
# tuning parameter
lambda = out$lambda
