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
  double max_pos = 0.0;  // Track maximum positive deviation
  double max_neg = 0.0;  // Track maximum negative deviation
  
  for(int j = 0; j < n; j++) {
    double step = isAbove[j] ? std::pow(scores[j], p_weight) / N_R : -1.0 / N_non_specific;
    running_sum += step;
    
    if(running_sum > max_pos) {
      max_pos = running_sum;
    }
    if(running_sum < max_neg) {
      max_neg = running_sum;
    }
  }
  
  // Return the max deviation with the larger absolute value
  return std::fabs(max_pos) > std::fabs(max_neg) ? max_pos : max_neg;
}
