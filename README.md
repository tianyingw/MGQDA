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
set.seed(1)
Xtrain <- matrix(rnorm(p * n * G), n * G, p)

# find projection matrix $\Omega$ with selected tuning parameter $\alpha = 0.1$ and $\lambda = 0.5$:
V <- dLDA(xtrain, ytrain, lambda = 0.1)
sum(rowSums(V) != 0)

# generate test data
m <- 20
set.seed(3)
xtest <- matrix(rnorm(p * m), m, p)

# perform classification
ytest <- classifyV(xtrain, ytrain, xtest, V)
