// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector calcRunningSum(NumericVector scores, LogicalVector isAbove, 
                           double N_R, double N_non_specific, double p_weight) {
  int n = scores.size();
  NumericVector rs(n);
  
  // First position
  if(isAbove[0]) {
    rs[0] = pow(scores[0], p_weight) / N_R;
  } else {
    rs[0] = -1.0 / N_non_specific;
  }
  
  // Remaining positions
  for(int j = 1; j < n; j++) {
    double step = isAbove[j] ? pow(scores[j], p_weight) / N_R : -1.0 / N_non_specific;
    rs[j] = rs[j-1] + step;
  }
  
  return rs;
}

// [[Rcpp::export]]
double calcMaxDeviation(NumericVector scores, LogicalVector isAbove, 
                      double N_R, double N_non_specific, double p_weight) {
  int n = scores.size();
  double running_sum = 0.0;
  double max_deviation = 0.0;
  
  for(int j = 0; j < n; j++) {
    double step = isAbove[j] ? pow(scores[j], p_weight) / N_R : -1.0 / N_non_specific;
    running_sum += step;
    
    if(fabs(running_sum) > max_deviation) {
      max_deviation = fabs(running_sum);
    }
  }
  
  return running_sum > 0 ? max_deviation : -max_deviation;
}