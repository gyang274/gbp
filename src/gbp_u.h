#ifndef __BINPACK_BINPACK_U__
#define __BINPACK_BINPACK_U__


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// using namespace arma;


// gbp_u.cpp

arma::umat create_max_mode_tbl(arma::umat& m, arma::uvec& ulmt, arma::uvec& vlmt, const int mval, const int mlvl);

// gbp_u.cpp - unique

arma::mat unique_rows(const arma::mat& m);

arma::mat unique_cols(const arma::mat& m);

// gbp_u.cpp - sort and sort_index

template <typename T>
inline bool approx_equal_cpp(const T& lhs, const T& rhs, double tol = 0.00000001) {
  return arma::approx_equal(lhs, rhs, "absdiff", tol);
}

arma::uvec sort_index_via_cols_internal(const arma::mat& m, const arma::uvec& ulmt, const arma::uvec& vlmt);

arma::uvec sort_index_via_cols(const arma::mat& m, const arma::uvec& vlmt);

arma::mat sort_via_cols(const arma::mat& m, const arma::uvec& vlmt);

arma::uvec sort_index_via_rows_internal(const arma::mat& m, const arma::uvec& ulmt, const arma::uvec& vlmt);

arma::uvec sort_index_via_rows(const arma::mat& m, const arma::uvec& ulmt);

arma::mat sort_via_rows(const arma::mat& m, const arma::uvec& ulmt);


#endif // __BINPACK_BINPACK_U__
