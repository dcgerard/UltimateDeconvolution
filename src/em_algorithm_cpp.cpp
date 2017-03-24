#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
#include <progress.hpp>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppProgress)]]

// Declare needed functions --------------------------------------------------
NumericMatrix get_llike_mat_cpp(const NumericMatrix& x_mat, const NumericMatrix& s_mat,
                                const NumericMatrix& v_mat, const NumericVector& pi_vec);
double dmixlike_cpp(const NumericMatrix& x_mat, const NumericMatrix& s_mat,
                    const NumericMatrix& v_mat, const NumericVector& pi_vec,
                    bool return_log = false);

// End Declariations ---------------------------------------------------------

//' Fixed point iteration from the EM algorithm.
//'
//' Note that I am changing v_mat and pi_vec by reference, but also returning them in the list.
//'
//' @inheritParams dmixlike
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
void em_fix_cpp(const NumericMatrix& x_mat, const NumericMatrix& s_mat,
                NumericMatrix& v_mat, NumericVector& pi_vec) {

  int N = x_mat.nrow();
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
  pi_vec = Rcpp::colSums(w_mat);
  pi_vec = pi_vec / Rcpp::sum(pi_vec);

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
  v_mat = Rcpp::wrap((lincom_s % lincom_xs).t());

  return;
}


//' C++ version of EM algorithm.
//'
//' @inheritParams ultimate_deconvolution
//' @param plot_iter A logical. Should we plot updates (\code{TRUE}) or not (\code{FALSE})?
//'
//' @return A list of the following elements:
//'
//'     \code{pi_vec}: The final estimate of the mixing proportions.
//'
//'     \code{v_mat}: The final estimate of the square roots of the rank-1 covariance matrices.
//'
//'     \code{llike_vec}: The vector of log-likelihoods. Should be increasing.
//'
//'     \code{convergence}: A value of \code{0} indicates convergence. A value of \code{1} indicates that
//'         the limit \code{itermax} has been reached. A vlue of \code{2} indicates that the user
//'         interupted the optimization.
//'
//' @author David Gerard
//'
//' @export
//'
// [[Rcpp::export]]
List em_cpp(const NumericMatrix& x_mat, const NumericMatrix& s_mat,
            NumericMatrix& v_mat, NumericVector& pi_vec,
            int itermax = 500, double tol = 10 ^ -5,
            bool plot_iter = false) {

  // Starting log-likelihood
  double llike_current = dmixlike_cpp(x_mat, s_mat, v_mat, pi_vec, true);
  int convergence;

  int iterindex = 0;
  double err = tol + 1.0;
  double llike_old;
  std::vector<double> llike_vec;
  llike_vec.reserve(itermax); // set aside at least itermax space for llike_vec
  llike_vec.push_back(llike_current);

  double tol_like_increase = -1.0e-12;

  // Progress and plotting pre-processing --------------------------------
  Progress p(itermax, plot_iter);

  while (iterindex < itermax && err > tol) {

    // Check progress interupt -------------------------------------------
    p.increment(); // update progress
    if (Progress::check_abort()) {
      convergence = 2;
      return List::create(_["pi_vec"] = pi_vec, _["v_mat"] = v_mat,
                          _["llike_vec"] = llike_vec, _["convergence"] = convergence);
    }

    // Make updates ------------------------------------------------------
    llike_old = llike_current;
    em_fix_cpp(x_mat, s_mat, v_mat, pi_vec);
    llike_current = dmixlike_cpp(x_mat, s_mat, v_mat, pi_vec, true);

    // Check llike increases ---------------------------------------------
    if (llike_current - llike_old < tol_like_increase) {
      throw std::domain_error("Likelihood did not increase.");
    }

    err = std::abs(std::exp(llike_current - llike_old) - 1);
    llike_vec.push_back(llike_current);
    iterindex++;
  }

  if (iterindex == itermax) {
    convergence = 1;
  } else {
    convergence = 0;
  }


  return List::create(_["pi_vec"] = pi_vec, _["v_mat"] = v_mat, _["llike_vec"] = llike_vec,
                      _["convergence"] = convergence);
}





