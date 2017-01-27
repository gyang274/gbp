#include "gbp3d_ck.h"

// class gbp3d {
//
// public:
//
//   arma::vec p; arma::mat it; arma::vec bn; arma::uvec k; double o; bool ok;
//
//   // .constructor
//   gbp3d(
//     arma::vec p, arma::mat it, arma::vec bn, arma::uvec k, double o, bool ok
//   ) : p(p), it(it), bn(bn), k(k), o(o), ok(ok) {
//   }
//
// };

// //' gbp3d_checkr
// //' @description
// //'  auxilium of gbp3d_solver_dpp
// //' @details
// //'  check fit solution is valid:
// //'   no conflict between item and bin, and no conflict between each pair of item.
// //' @param sn <gbp3d>
// //'  gbp3d object from gbp3d_solver_dpp() solution.
// //' @return okfit? <bool>
// //' @family gbp3d_solver_dpp
// //' @export
// // [[Rcpp::export]]
bool gbp3d_checkr(gbp3d sn) {

  bool okfit = true;

  arma::uvec klmt = arma::find(sn.k == 1);

  // main
  // main: no conflict between it and bn
  for (arma::uword i = 0; i < klmt.size(); i++) {
    if (sn.it(0, klmt(i)) + sn.it(3, klmt(i)) > sn.bn(0) ||
        sn.it(1, klmt(i)) + sn.it(4, klmt(i)) > sn.bn(1) ||
        sn.it(2, klmt(i)) + sn.it(5, klmt(i)) > sn.bn(2)
    ) {
      Rcpp::Rcout << "gbp3d_checkr: it conflict bn: index " << klmt(i) << " ." << std::endl;
      okfit = false; return okfit;
    }
  }

  // main: no conflict between each pair of it
  for (arma::uword i = 0; i < klmt.size(); i++) {
    for (arma::uword j = i + 1; j < klmt.size(); j++) {
      if (!(
          sn.it(0, klmt(i)) + sn.it(3, klmt(i)) <= sn.it(0, klmt(j)) || sn.it(0, klmt(j)) + sn.it(3, klmt(j)) <= sn.it(0, klmt(i)) ||
          sn.it(1, klmt(i)) + sn.it(4, klmt(i)) <= sn.it(1, klmt(j)) || sn.it(1, klmt(j)) + sn.it(4, klmt(j)) <= sn.it(1, klmt(i)) ||
          sn.it(2, klmt(i)) + sn.it(5, klmt(i)) <= sn.it(2, klmt(j)) || sn.it(2, klmt(j)) + sn.it(5, klmt(j)) <= sn.it(2, klmt(i))
      )) {
        Rcpp::Rcout << "gbp3d_checkr: it conflict it: index " << klmt(i) << " and " << klmt(j) << "." << std::endl;
        okfit = false; return okfit;
      }
    }
  }

  return okfit;
}


// class gbp3q {
//
// public:
//
//   arma::vec p; arma::mat it; arma::mat bn; arma::uvec k; arma::uvec f; double o; bool ok;
//
//   // .constructor
//   gbp3q(
//     arma::vec p, arma::mat it, arma::mat bn, arma::uvec k, arma::uvec f, double o, bool ok
//   ) : p(p), it(it), bn(bn), k(k), f(f), o(o), ok(ok) {
//   }
//
// };

// //' gbp3q_checkr
// //' @description
// //'  auxilium of gbp3d_solver_dpp_filt
// //' @details
// //'  check fit solution is valid:
// //'   no conflict between item and bin, and no conflict between each pair of item.
// //' @param sn <gbp3q>
// //'  gbp3q object from gbp3d_solver_dpp_filt() solution.
// //' @return okfit? <bool>
// //' @family gbp3d_solver_dpp_filt
// //' @export
// // [[Rcpp::export]]
bool gbp3q_checkr(gbp3q sn) {

  bool okfit = false;

  arma::uvec id = arma::find(sn.f == 1);

  if (id.size() == 1) {
    okfit = gbp3d_checkr(gbp3d(sn.p, sn.it, sn.bn.col(id(0)), sn.k, sn.o, sn.ok));
  } else {
    Rcpp::Rcout << "gbp3q_checkr: f should have a unique index label 1." << std::endl;
  }

  return okfit;
}

