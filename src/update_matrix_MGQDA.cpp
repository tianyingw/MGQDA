// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// 定义 f 函数  
// [[Rcpp::export]]
double f(double x, double a, double b, double alpha, double lambda) {  
  double sum = a * x + alpha * lambda * x / std::sqrt(std::pow(b, 2) + std::pow(x, 2));  
  return sum;  
}  

// 定义 BiSec 函数  
// [[Rcpp::export]]
double BiSec(double mi, double ma, double a, double b, double c, double alpha, double lambda) {  
  double epsilon = 1e-8;  
  while (true) {  
    double xm = (ma + mi) / 2;  
    if (f(xm, a, b, alpha, lambda) < c) {  
      mi = xm;  
    } else {  
      ma = xm;  
    }  
    if (abs(ma - mi) < epsilon) {  
      break;  
    }  
  }  
  return ma;  
}  


// 定义一个函数，用于更新矩阵V并返回结果
// S, B, D是传入的矩阵参数
// G, p, alpha, lambda是其他参数
// eps, maxiter是收敛条件
// 返回一个列表，包含更新后的V矩阵和nfeature值
// [[Rcpp::depends(RcppArmadillo)]]
List updateMatrix(mat S, mat B, mat D, mat V,
                  double lambda, double alpha, double eps, int maxiter, int p, int G) {
  int N = 0;
  
  // 定义中间变量
  mat V0(p, G*(G-1));
  rowvec Vec(G-1);
  double a, b, c, x;
  
  // 迭代更新矩阵
  do {
    V0 = V;
    for (int j = 0; j < p; ++j) {
      for (int k = 0; k < G; ++k) {
        Vec = (S(j, span(p*k, p*(k+1)-1)) + B.row(j)) * V.cols(span((G-1)*k, (G-1)*(k+1)-1)) - D.row(j);
        Vec = Vec - (S(j, p*k+j) + B(j,j)) * V(j, span((G-1)*k, (G-1)*(k+1)-1));
        
        a = S(j, p*k+j) + B(j,j);
        b = std::sqrt(sum(V.row(j) % V.row(j)) - sum(V(j, span((G-1)*k, (G-1)*(k+1)-1)) % V(j, span((G-1)*k, (G-1)*(k+1)-1))));
        c = std::sqrt(sum(Vec % Vec));
        
        if (c < (lambda*(1-alpha)/(std::sqrt(G)))) {
          V(j, span((G-1)*k, (G-1)*(k+1)-1)) = rowvec((G-1), fill::zeros);  // 清零
        } else {
          c = c - lambda*(1-alpha)/(std::sqrt(G));
          x = BiSec(0.0, c/a, a, b, c, alpha, lambda);
          V(j, span((G-1)*k, (G-1)*(k+1)-1)) = -x*Vec/c;
        }
      }
    }
    N++;
  } while (arma::abs(V0-V).max() >= eps && N <= maxiter);
  
  // 处理 V 矩阵中的较小值，将其设为零
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < (G*(G-1)); ++j) {
      if (abs(V(i, j)) < 1e-6) {
        V(i, j) = 0.0;
      }
    }
  }
  // 返回结果列表
  return List::create(Named("V") = V,  Named("iter") = N);
}


// [[Rcpp::depends(RcppArmadillo)]]
List updateMatrix_2(mat S, mat B, vec D, mat V,
                  double lambda, double alpha, double eps, int maxiter, int p) {
  int N = 0;
  
  // 定义中间变量
  mat V0(p, 2);
  double Vec, a, b, c, x;
  
  // 迭代更新矩阵
  do {
    V0 = V;
    for (int j = 0; j < p; ++j) {
      for (int k = 0; k < 2; ++k) {
        Vec = sum((S(j, span(p*k, p*(k+1)-1)) + B.row(j)) % V.col(k)) - D(j);
        Vec = Vec - (S(j, p*k+j) + B(j,j)) * V(j, k);
        
        a = S(j, p*k+j) + B(j,j);
        b = std::sqrt(sum(V.row(j) % V.row(j)) - (V(j, k) * V(j, k)));
        c = std::sqrt(Vec * Vec);
        
        if (c < (lambda*(1-alpha)/(std::sqrt(2)))) {
          V(j, k) = 0.0;  // 清零
        } else {
          c = c - lambda*(1-alpha)/(std::sqrt(2));
          x = BiSec(0.0, c/a, a, b, c, alpha, lambda);
          V(j, k) = -x*Vec/c;
        }
      }
    }
    N++;
  } while (arma::abs(V0-V).max() >= eps && N <= maxiter);
  
  // 处理 V 矩阵中的较小值，将其设为零
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < 2; ++j) {
      if (abs(V(i, j)) < 1e-6) {
        V(i, j) = 0.0;
      }
    }
  }
  // 返回结果列表
  return List::create(Named("V") = V,  Named("iter") = N);
}


// 导出C++函数，使其可以在R中调用
RCPP_MODULE(matrix_mult_module) {
  function("f", &f);
  function("BiSec", &BiSec);
  function("updateMatrix", &updateMatrix);
  function("updateMatrix_2", &updateMatrix_2);
}
