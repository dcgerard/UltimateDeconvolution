#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' Multivariate normal density at mean 0 and a covariance that is the sum
//' of a rank-1 matrix and a diagonal matrix.
//'
//' @inheritParams dr1_norm
//' @param return_log Should we return the log-density (\code{TRUE}) or not (\code{FALSE})?
//'
//' @author David Gerard
//'
//' @export
//'
// [[Rcpp::export]]
double dnorm_rank1(NumericVector x, NumericVector v, NumericVector s_diag,
                   NumericVector mu, bool return_log = false) {
  int R = x.size();
  double vsv = Rcpp::sum(v * v / s_diag);
  double xmusv = Rcpp::sum((x - mu) * v / s_diag);
  double xmusxmu = Rcpp::sum(Rcpp::pow(x - mu, 2) / s_diag);

  double llike = (-R * std::log(2 * PI) - std::log(1 + vsv) - Rcpp::sum(Rcpp::log(s_diag)) -
                  xmusxmu + std::pow(xmusv, 2) / (1 + vsv)) / 2;

  if (return_log) {
    return llike;
  }
  else {
    return std::exp(llike);
  }
}

