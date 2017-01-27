#include "gbp2d_ck.h"

// class gbp2d {
//
// public:
//
//   arma::vec p; arma::mat it; arma::vec bn; arma::uvec k; double o; bool ok;
//
//   // .constructor
//   gbp2d(
//     arma::vec p, arma::mat it, arma::vec bn, arma::uvec k, double o, bool ok
//   ) : p(p), it(it), bn(bn), k(k), o(o), ok(ok) {
//   }
//
// };

// //' gbp2d_checkr
// //' @description
// //'  auxilium of gbp2d_solver_dpp
// //' @details
// //'  check fit solution is valid:
// //'   no conflict between item and bin, and no conflict between each pair of item.
// //' @param sn <gbp2d>
// //'  gbp2d object from gbp2d_solver_dpp() solution.
// //' @return okfit? <bool>
// //' @family gbp2d_solver_dpp
// //' @export
// // [[Rcpp::export]]
bool gbp2d_checkr(gbp2d sn) {

  bool okfit = true;

  arma::uvec klmt = arma::find(sn.k == 1);

  // main
  // main: no conflict between it and bn
  for (arma::uword i = 0; i < klmt.size(); i++) {
    if (sn.it(0, klmt(i)) + sn.it(2, klmt(i)) > sn.bn(0) ||
        sn.it(1, klmt(i)) + sn.it(3, klmt(i)) > sn.bn(1)
    ) {
      Rcpp::Rcout << "gbp2d_checkr: it conflict bn: index " << klmt(i) << " ." << std::endl;
      okfit = false; return okfit;
    }
  }

  // main: no conflict between each pair of it
  for (arma::uword i = 0; i < klmt.size(); i++) {
    for (arma::uword j = i + 1; j < klmt.size(); j++) {
      if (!(
          sn.it(0, klmt(i)) + sn.it(2, klmt(i)) <= sn.it(0, klmt(j)) || sn.it(0, klmt(j)) + sn.it(2, klmt(j)) <= sn.it(0, klmt(i)) ||
          sn.it(1, klmt(i)) + sn.it(3, klmt(i)) <= sn.it(1, klmt(j)) || sn.it(1, klmt(j)) + sn.it(3, klmt(j)) <= sn.it(1, klmt(i))
      )) {
        Rcpp::Rcout << "gbp2d_checkr: it conflict it: index " << klmt(i) << " and " << klmt(j) << "." << std::endl;
        okfit = false; return okfit;
      }
    }
  }

  return okfit;
}


// class gbp2q {
//
// public:
//
//   arma::vec p; arma::mat it; arma::mat bn; arma::uvec k; arma::uvec f; double o; bool ok;
//
//   // .constructor
//   gbp2q(
//     arma::vec p, arma::mat it, arma::mat bn, arma::uvec k, arma::uvec f, double o, bool ok
//   ) : p(p), it(it), bn(bn), k(k), f(f), o(o), ok(ok) {
//   }
//
// };

// //' gbp2q_checkr
// //' @description
// //'  auxilium of gbp2d_solver_dpp_filt
// //' @details
// //'  check fit solution is valid:
// //'   no conflict between item and bin, and no conflict between each pair of item.
// //' @param sn <gbp2q>
// //'  gbp2q object from gbp2d_solver_dpp_filt() solution.
// //' @return okfit? <bool>
// //' @family gbp2d_solver_dpp_filt
// //' @export
// // [[Rcpp::export]]
bool gbp2q_checkr(gbp2q sn) {

  bool okfit = false;

  arma::uvec id = arma::find(sn.f == 1);

  if (id.size() == 1) {
    okfit = gbp2d_checkr(gbp2d(sn.p, sn.it, sn.bn.col(id(0)), sn.k, sn.o, sn.ok));
  } else {
    Rcpp::Rcout << "gbp2q_checkr: f should have a unique index label 1." << std::endl;
  }

  return okfit;
}

