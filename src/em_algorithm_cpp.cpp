#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// Declare needed functions
NumericMatrix get_llike_mat_cpp(const NumericMatrix& x_mat, const NumericMatrix& s_mat,
                                const NumericMatrix& v_mat, const NumericVector& pi_vec);

//' Fixed point iteration from the EM algorithm.
//'
//' @inheritParams dmixlike
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
List em_fix_cpp(const NumericMatrix& x_mat, const NumericMatrix& s_mat,
                const NumericMatrix& v_mat, const NumericVector& pi_vec) {

  int N = x_mat.nrow();
  int R = x_mat.ncol();
  int K = v_mat.ncol();

  // Get W matrix (the element-specific mixing proportions)
  NumericMatrix llike_mat = get_llike_mat_cpp(x_mat, s_mat, v_mat, pi_vec);

  NumericMatrix w_mat(N, K);
  double max_element;
  double ldenom;

  for (int obs_index = 0; obs_index < N; obs_index++) {
    max_element = Rcpp::max(llike_mat(obs_index, _));
    ldenom = std::log(Rcpp::sum(Rcpp::exp(llike_mat(obs_index, _) - max_element))) + max_element;
    w_mat(obs_index, _) = Rcpp::exp(llike_mat(obs_index, _) - ldenom);
  }

  // Update the pi-values ----------------------------------------
  NumericVector pi_new = Rcpp::colSums(w_mat);
  pi_new = pi_new / Rcpp::sum(pi_new);

  Rcpp::Rcout << pi_new << std::endl;

  // Get theta and eta matrices
  arma::mat theta_mat(N, K);
  arma::mat eta_mat(N, K);
  double sigma2_kj;
  double mu_kj;
  double vsv;
  double xsv;

  for (int obs_index = 0; obs_index < N; obs_index++) {
    for (int mix_index = 0; mix_index < K; mix_index++) {
      vsv = Rcpp::sum(Rcpp::pow(v_mat(_, mix_index), 2) / s_mat(obs_index, _));
      xsv = Rcpp::sum(v_mat(_, mix_index) * x_mat(obs_index, _) / s_mat(obs_index, _));

      sigma2_kj = 1 / (vsv + 1);
      mu_kj = sigma2_kj * xsv;

      theta_mat(obs_index, mix_index) = w_mat(obs_index, mix_index) * mu_kj;
      eta_mat(obs_index, mix_index) = w_mat(obs_index, mix_index) *
        (std::pow(mu_kj, 2) + sigma2_kj);
    }
  }

  // Update v_mat
  arma::mat s_arma = as<arma::mat>(s_mat); // NB: s_arma points to same space in memory as s_mat
  arma::mat x_arma = as<arma::mat>(x_mat); // NB: x_arma points to same space in memory as x_mat
  arma::mat lincom_s = 1 / (eta_mat.t() * (1 / s_arma));
  arma::mat lincom_xs = theta_mat.t() * (x_arma / s_arma);
  arma::mat v_new = (lincom_s % lincom_xs).t();

  return List::create(_["pi_vec"] = pi_new, _["v_mat"] = v_new);
}







