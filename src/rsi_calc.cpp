#include <Rcpp.h>
#include <vector>        // for std::vector
#include <functional>    // for std::less, etc
#include <algorithm>     // for count_if

using namespace Rcpp ;

// [[Rcpp::export]]
int rsi_calc_S(std::vector<double> x, bool include_I) {
  if (include_I == TRUE) {
    return count_if(x.begin(), x.end(), bind2nd(std::less_equal<double>(), 2));
  } else {
    return count_if(x.begin(), x.end(), bind2nd(std::less<double>(), 2));
  }
}

// [[Rcpp::export]]
int rsi_calc_R(std::vector<double> x, bool include_I) {
  if (include_I == TRUE) {
    return count_if(x.begin(), x.end(), bind2nd(std::greater_equal<double>(), 2));
  } else {
    return count_if(x.begin(), x.end(), bind2nd(std::greater<double>(), 2));
  }
}

// [[Rcpp::export]]
int rsi_calc_total(std::vector<double> x) {
 return count_if(x.begin(), x.end(), bind2nd(std::less_equal<double>(), 3));
}
