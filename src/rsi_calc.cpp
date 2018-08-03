#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
int rsi_calc_S(DoubleVector x, bool include_I) {
  return count_if(x.begin(),
                  x.end(),
                  bind2nd(std::less_equal<double>(),
                          1 + include_I));
}

// [[Rcpp::export]]
int rsi_calc_I(DoubleVector x) {
  return count_if(x.begin(),
                  x.end(),
                  bind2nd(std::equal_to<double>(),
                          2));
}

// [[Rcpp::export]]
int rsi_calc_R(DoubleVector x, bool include_I) {
  return count_if(x.begin(),
                  x.end(),
                  bind2nd(std::greater_equal<double>(),
                          3 - include_I));
}
