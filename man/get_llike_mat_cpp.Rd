% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{get_llike_mat_cpp}
\alias{get_llike_mat_cpp}
\title{Get's a matrix of likeklihood values.}
\usage{
get_llike_mat_cpp(x_mat, s_mat, v_mat, pi_vec)
}
\arguments{
\item{x_mat}{A matrix of data. The rows index the observations and
the columns index the variables.}

\item{s_mat}{A matrix of variances (NOT standard deviations).
The rows index the observations and the columns index the
variables.}

\item{v_mat}{A matrix of square-roots of the rank-1 mixing
covariance matrices. The columns index the mixing
components, the rows index the variables.}

\item{pi_vec}{A vector of mixture proportions.}
}
\description{
Element (i, j) of the returned matrix is the log of
pi_j N(x_i | 0, v_j v_j^T + S_i)
}
\author{
David Gerard
}
