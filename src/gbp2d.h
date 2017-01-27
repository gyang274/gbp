#ifndef __BINPACK_BINPACK2D__
#define __BINPACK_BINPACK2D__


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// using namespace arma;


// gbp2d
// @description
//  generalized bin packing problem in 2 dimension, a.k.a rectangle fill.
// @details
//  gbp2d init a profit vector p, a length vector l, a depth vector d,
//   a length constraint ml, and a depth constraint md on l x d rectangle
//   with geometry intepretation.
//
//  gbp2d solver would solve
//
//    maximize   sum_{j=1}^{n} p_{j} k_{j}
//
//    subject to fit (l_{j}, d_{j}) at coordinate (x_{j}, y_{j})
//               such that no overlap in ml x md, j = 1, ...., n
//
//  and instantiate a gbp2d object with a x-axis coordinate vector x, a y-axis coordinate vector y,
//   a selection vector k, and an objective o.
// @note
//  p is a proxy of ranking on rectangle fit difficulty, often a function w.r.t max(l, d) and l x d
//   and solver would often maximize sum_{j=1}^{n} v_{j} k_{j} instead of sum_{j=1}^{n} p_{j} k_{j}
// @family gbp2d
// @rdname gbp2d
// class gbp2d {
//
// public:
//
//   arma::vec p;
//
//   arma::uvec l; arma::uvec d; arma::uword ml; arma::uword md;
//
//   arma::uvec x; arma::uvec y; arma::uvec k; double o;
//
//   // .constructor
//   gbp2d(
//     arma::vec p, arma::umat ld, arma::uvec m, arma::umat xy, arma::uvec k, double o
//   ) : p(p), k(k), o(o) {
//     l = ld.row(0); d = ld.row(1); ml = m(0); md = m(1);
//     x = xy.row(0); y = xy.row(1);
//   }
//
// };
//
// gbp2d
// @param p
//  p profit of it fit into bn <vector>
//  - cluster max(l, d) and min(l, d) via gbp2d_solver_dpp_prep_create_p()
// @param it
//  it position and scale <matrix>
//  - x, y it position in the bin <numeric>
//  - l, d it scale along x and y <numeric>
// @param bn
//  bn scale <vector>
//  - l, d bn scale along x and y <numeric>
// @param k
//  k selection indicator 0, 1 <vector>
// @param o
//  o objective achivement volumn fit in over volumn overall <numeric>
// @param ok
//  ok a quick indicator of all it fit into bn? <bool>
class gbp2d {

public:

  arma::vec p; arma::mat it; arma::vec bn; arma::uvec k; double o; bool ok;

  // .constructor
  gbp2d(
    arma::vec p, arma::mat it, arma::vec bn, arma::uvec k, double o, bool ok
  ) : p(p), it(it), bn(bn), k(k), o(o), ok(ok) {
  }

};


// solve gbp2d via extreme point heuristic and best information score fit strategy
gbp2d gbp2d_solver_dpp(const arma::vec& p, const arma::mat& ld, const arma::vec& m);

bool gbp2d_solver_dpp_main(
  const arma::vec& bn, arma::mat& it, arma::mat& itastr, arma::mat& xp, const arma::uvec& q, const arma::uword nlvl, const arma::uword nastr,
  arma::uvec g, arma::uvec& gastr, const arma::vec& v, const double& vastr, double u, double& uastr
);

arma::uword gbp2d_solver_dpp_main_create_nastr(const arma::vec& p, const arma::mat& ld, const arma::vec& m);

arma::uword gbp2d_solver_dpp_main_create_nlmt(const arma::uword nlvl, const arma::uword nastr);

arma::vec gbp2d_solver_dpp_prep_create_p(const arma::mat& ld, const arma::vec& m);


// gbp2q
// @param p
//  p profit of it fit into bn <vector>
//  - cluster max(l, d) and min(l, d) via gbp2d_solver_dpp_main_create_p()
// @param it
//  it position and scale <matrix>
//  - x, y it position in the bin <numeric>
//  - l, d it scale along x and y <numeric>
// @param bn
//  bn scale <matrix>
//  - l, d bn scale along x and y <numeric>
//  - l, d in row and each col is a single bn
//  often the first col is the most prefered smallest bn while should
//   the last col is the least prefered largest and often dominant bn
//  should make sure no X in front of Y if bnX dominant bnY,
//   bnX dominant bnY if all(X(l, d) > Y(l, d)) and should always prefer Y.
//  should make sure bn such that l >= d or vice versa.
// @param k
//  k selection indicator 0, 1 on it <vector>
// @param f
//  f selection indicator 0, 1, 2, 3 on bn <vector>
//  f in result should have no 0 and only one of 1
// @param o
//  o objective achivement volumn fit in over volumn overall <numeric>
// @param ok
//  ok a quick indicator of all it fit into bn? <bool>
class gbp2q {

public:

  arma::vec p; arma::mat it; arma::mat bn; arma::uvec k; arma::uvec f; double o; bool ok;

  // .constructor
  gbp2q(
    arma::vec p, arma::mat it, arma::mat bn, arma::uvec k, arma::uvec f, double o, bool ok
  ) : p(p), it(it), bn(bn), k(k), f(f), o(o), ok(ok) {
  }

};


// solve gbp2d w.r.t select most preferable often smallest bin from bn list
gbp2q gbp2d_solver_dpp_filt(const arma::mat& ld, const arma::mat& m);

// bool assert_matrix_last_col_dominant(const arma::mat& m);

void gbp2d_solver_dpp_filt_fast(
  const arma::mat& sld, const double vld, const arma::mat& sm, const arma::rowvec& vm, arma::uvec& flmt
);

void gbp2d_solver_dpp_filt_slow(
  const double vld, const arma::rowvec& vm, arma::uvec& flmt
);

void gbp2d_solver_dpp_filt_dcol(
  const arma::mat &sm, const arma::uword icol, arma::uvec& flmt
);


#endif // __BINPACK_BINPACK2D__
