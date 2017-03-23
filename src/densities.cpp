#include <Rcpp.h>
#include <math.h>
#include <cassert>
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

//' Get's a matrix of likeklihood values.
//'
//' Element (i, j) of the returned matrix is the log of
//' pi_j N(x_i | 0, v_j v_j^T + S_i)
//'
//' @inheritParams dmixlike
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
NumericMatrix get_llike_mat_cpp(NumericMatrix x_mat, NumericMatrix s_mat,
                                NumericMatrix v_mat, NumericVector pi_vec) {

  // Check input -------------------------------------------------------------
  int R = x_mat.ncol();
  int N = x_mat.nrow();
  int K = v_mat.ncol();

  assert(R == s_mat.ncol());
  assert(N == s_mat.nrow());
  assert(K == pi_vec.siz());
  assert(R == v_mat.nrow());
  assert(std::abs(Rcpp::sum(pi_vec) - 1) < 10 ^ -10);

  NumericMatrix llike_mat(N, K);

  NumericVector mu(R); // All zeros

  double llike;

  NumericVector lpi = Rcpp::log(pi_vec);

  for (int obs_index = 0; obs_index < N; obs_index++) {
    for (int mix_index = 0; mix_index < K; mix_index++) {
      llike = dnorm_rank1(x_mat(obs_index, _), v_mat(_, mix_index),
                          s_mat(obs_index, _), mu, true);
      llike_mat(obs_index, mix_index) = llike + lpi(mix_index);
    }
  }

  return(llike_mat);
}

