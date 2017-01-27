#ifndef __BINPACK_BINPACK3D__
#define __BINPACK_BINPACK3D__


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// using namespace arma;


// gbp3d
// @description
//  generalized bin packing problem in 3 dimension, a.k.a bin packing problem.
// @details
//  gbp3d init a profit vector p, a length vector l, a depth vector d, a height vector h, and also
//   a length constraint ml, a depth constraint md, and a height constraint mh on l x d x h cuboid
//   with geometry intepretation.
//
//  gbp3d solver would solve
//
//    maximize   sum_{j=1}^{n} p_{j} k_{j}
//
//    subject to fit (l_{j}, d_{j}, h_{j}) at coordinate (x_{j}, y_{j}, z_{j})
//               such that no overlap in ml x md x mh cuboid, j = 1, ......, n
//
//  and instantiate a gbp3d object with a x-axis coordinate vector x, a y-axis coordinate vector y,
//   a z-axis coordinate vector z, a selection vector k, and an objective o.
// @note
//  p is a proxy of ranking on cuboid fit difficulty, often a func of max(l, d, h), surface, volume
//   and solver would often maximize sum_{j=1}^{n} v_{j} k_{j} instead of sum_{j=1}^{n} p_{j} k_{j}
// @family gbp3d
// @rdname gbp3d
// class gbp3d {
//
// public:
//
//   arma::vec p;
//
//   arma::vec l; arma::vec d; arma::vec h; double ml; double md; double mh;
//
//   arma::vec x; arma::vec y; arma::vec z; arma::uvec k; double o;
//
//   // .constructor
//   gbp3d(
//     arma::vec p, arma::mat ldh, arma::vec m, arma::mat xyz, arma::uvec k, double o
//   ) : p(p), k(k), o(o) {
//     l = ldh.row(0); d = ldh.row(1); h = ldh.row(2); ml = m(0); md = m(1); mh = m(2);
//     x = xyz.row(0); y = xyz.row(1); z = xyz.row(2);
//   }
//
// };
//
// gbp3d
// @param p
//  p profit of it fit into bn <vector>
//  - cluster max(l,d,h) and area via gbp3d_solver_dpp_main_create_p()
// @param it
//  it position and scale <matrix>
//  - x, y, z it position in the bin <numeric>
//  - l, d, h it scale along x, y, z <numeric>
// @param bn
//  bn scale <vector>
//  - l, d, h bn scale along x, y, z <numeric>
// @param k
//  k selection indicator 0, 1 <vector>
// @param o
//  o objective achivement volumn fit in over volumn overall <numeric>
// @param ok
//  ok a quick indicator of all it fit into bn? <bool>
class gbp3d {

public:

  arma::vec p; arma::mat it; arma::vec bn; arma::uvec k; double o; bool ok;

  // .constructor
  gbp3d(
    arma::vec p, arma::mat it, arma::vec bn, arma::uvec k, double o, bool ok
  ) : p(p), it(it), bn(bn), k(k), o(o), ok(ok) {
  }

};


// solve gbp3d via extreme point heuristic and best information score fit strategy
gbp3d gbp3d_solver_dpp(const arma::vec& p, const arma::mat& ldh, const arma::vec& m);

bool gbp3d_solver_dpp_main(
  const arma::vec& bn, arma::mat& it, arma::mat& itastr, arma::mat& xp, const arma::uvec& q, const arma::uword nlvl, const arma::uword nastr,
  arma::uvec g, arma::uvec& gastr, const arma::vec& v, const double& vastr, double u, double& uastr
);

arma::uword gbp3d_solver_dpp_main_create_nastr(const arma::vec& p, const arma::mat& ldh, const arma::vec& m);

arma::uword gbp3d_solver_dpp_main_create_nlmt(const arma::uword nlvl, const arma::uword nastr);

arma::vec gbp3d_solver_dpp_prep_create_p(const arma::mat& ldh, const arma::vec& m);


// gbp3q
// @param p
//  p profit of it fit into bn <vector>
//  - cluster max(l,d,h) and area via gbp3d_solver_dpp_main_create_p()
// @param it
//  it position and scale <matrix>
//  - x, y, z it position in the bin <numeric>
//  - l, d, h it scale along x, y, z <numeric>
// @param bn
//  bn scale <matrix>
//  - l, d, h bn scale along x, y, z <numeric>
//  - l, d, h in row and each col is a single bn
//  should have the first column with the most prefered smallest bn while having
//   the last column as the least prefered largest bn
//  should make sure no X in front of Y if bnX dominant bnY,
//   bnX dominant bnY if all(X(l, d, h) > Y(l, d, h)) and should always prefer Y.
//  should make sure bn such that l >= d >= h or vice versa.
// @param k
//  k selection indicator 0, 1 on it <vector>
// @param f
//  f selection indicator 0, 1, 2, 3 on bn <vector>
//  f in result should have no 0 and only one of 1
// @param o
//  o objective achivement volumn fit in over volumn overall <numeric>
// @param ok
//  ok a quick indicator of all it fit into bn? <bool>
class gbp3q {

public:

  arma::vec p; arma::mat it; arma::mat bn; arma::uvec k; arma::uvec f; double o; bool ok;

  // .constructor
  gbp3q(
    arma::vec p, arma::mat it, arma::mat bn, arma::uvec k, arma::uvec f, double o, bool ok
  ) : p(p), it(it), bn(bn), k(k), f(f), o(o), ok(ok) {
  }

};


// solve gbp3d w.r.t select most preferable often smallest bin from bn list
gbp3q gbp3d_solver_dpp_filt(const arma::mat& ldh, const arma::mat& m);

// bool assert_matrix_last_col_dominant(const arma::mat& m);

void gbp3d_solver_dpp_filt_fast(
  const arma::mat& sldh, const double vldh, const arma::mat& sm, const arma::rowvec& vm, arma::uvec& flmt
);

void gbp3d_solver_dpp_filt_slow(
  const double vldh, const arma::rowvec& vm, arma::uvec& flmt
);

void gbp3d_solver_dpp_filt_dcol(
  const arma::mat &sm, const arma::uword icol, arma::uvec& flmt
);


#endif // __BINPACK_BINPACK3D__
