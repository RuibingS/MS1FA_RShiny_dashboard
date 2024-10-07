#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Function to calculate correlation based on the specified method
double calculateCorrelation(NumericVector col1, NumericVector col2, std::string method) {
  LogicalVector completeIndices = !is_na(col1) & !is_na(col2);
  
  if (sum(completeIndices) > 1) {
    arma::vec v1 = as<arma::vec>(col1[completeIndices]);
    arma::vec v2 = as<arma::vec>(col2[completeIndices]);
    
    if (method == "pearson") {
      return arma::as_scalar(arma::cor(v1, v2));
    } else if (method == "spearman" || method == "kendall") {
      Function order("order");
      NumericVector rank1 = wrap(order(v1));
      NumericVector rank2 = wrap(order(v2));
      
      Function cor("cor");
      return as<double>(cor(rank1, rank2, Named("method") = method));
    } else {
      stop("Unsupported correlation method. Use 'pearson', 'spearman', or 'kendall'.");
    }
  } else {
    return NA_REAL;
  }
}

// [[Rcpp::export]]
NumericMatrix pairwiseCor(NumericMatrix x, std::string method) {
  int nCols = x.ncol();
  NumericMatrix result(nCols, nCols);
  
  for (int i = 0; i < nCols - 1; ++i) {
    for (int j = i + 1; j < nCols; ++j) {
      NumericVector col1 = x(_, i);
      NumericVector col2 = x(_, j);
      
      // Calculate correlation based on the specified method
      result(i, j) = calculateCorrelation(col1, col2, method);
      result(j, i) = result(i, j);  // Fill the lower triangular part of the matrix
    }
  }
  
  // Set diagonal elements to 1
  for (int i = 0; i < nCols; ++i) {
    result(i, i) = 1.0;
  }
  
  return result;
}
