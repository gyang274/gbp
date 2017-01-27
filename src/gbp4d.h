#ifndef __BINPACK_BINPACK4D__
#define __BINPACK_BINPACK4D__


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// using namespace arma;


// gbp4d
// @description
//  generalized bin packing problem in 4 dimension, a.k.a bin packing problem with weight limit.
// @details
//  gbp4d init a profit vector p, a length l, a depth d, a height h, and a weight w, along with
//   associate constraints ml, md, mh and mw.
//  gbp4d should fit it (l, d, h, w) into bn (ml, md, mh, mw) with w on weight limit constraint
//   and l, d, h on geometry intepretation.
//  gbp4d solver would solve
//
//    maximize   sum_{j=1}^{n} p_{j} k_{j}
//
//    subject to sum_{j=1}^{n} w_{j} k_{j} leq mw and
//
//               fit (l_{j}, d_{j}, h_{j}) at coordinate (x_{j}, y_{j}, z_{j})
//               such that no overlap in ml x md x mh cuboid, j = 1, ......, n
//
//  and instantiate a gbp4d object with a x-axis coordinate vector x, a y-axis coordinate vector y,
//   a z-axis coordinate vector z, a selection vector k, and an objective o.
// @note
//  p is a proxy of ranking on cuboid fit difficulty, often a func of max(l, d, h), surface, volume
//   and solver would often maximize sum_{j=1}^{n} v_{j} k_{j} instead of sum_{j=1}^{n} p_{j} k_{j}
// @family gbp4d
// @rdname gbp4d
// class gbp4d {
//
// public:
//
//   arma::vec p;
//
//   arma::vec l; arma::vec d; arma::vec h; arma::vec w; double ml; double md; double mh; double mw;
//
//   arma::vec x; arma::vec y; arma::vec z; arma::vec w; arma::uvec k; double o;
//
//   // .constructor
//   gbp3d(
//     arma::vec p, arma::mat ldhw, arma::vec m, arma::mat xyzw, arma::uvec k, double o
//   ) : p(p), k(k), o(o) {
//     l = ldhw.row(0); d = ldhw.row(1); h = ldhw.row(2); w = ldhw.row(3); ml = m(0); md = m(1); mh = m(2); mw = m(3);
//     x = xyzw.row(0); y = xyzw.row(1); z = xyzw.row(2); w = xyzw.row(3);
//   }
//
// };
//
// gbp4d
// @param p
//  p profit of it fit into bn <vector>
//  - cluster w, cluster max(l,d,h) and area via gbp4d_solver_dpp_main_create_p()
// @param it
//  it position and scale <matrix>
//  - x, y, z, w it position and w hold in the bin <numeric>
//  - l, d, h, w it scale along x, y, z and also w <numeric>
// @param bn
//  bn scale <vector>
//  - l, d, h, w bn scale along x, y, z and also w <numeric>
// @param k
//  k selection indicator 0, 1 <vector>
// @param o
//  o objective achivement volumn fit in over volumn overall <numeric>
// @param ok
//  ok a quick indicator of all it fit into bn? <bool>
class gbp4d {

public:

  arma::vec p; arma::mat it; arma::vec bn; arma::uvec k; double o; bool ok;

  // .constructor
  gbp4d(
    arma::vec p, arma::mat it, arma::vec bn, arma::uvec k, double o, bool ok
  ) : p(p), it(it), bn(bn), k(k), o(o), ok(ok) {
  }

};


// solve gbp4d via extreme point heuristic and best information score fit strategy
gbp4d gbp4d_solver_dpp(const arma::vec& p, const arma::mat& ldhw, const arma::vec& m);

bool gbp4d_solver_dpp_main(
    const arma::vec& bn, arma::mat& it, arma::mat& itastr, arma::mat& xp, const arma::uvec& q, const arma::uword nlvl, const arma::uword nastr,
    arma::uvec g, arma::uvec& gastr, const arma::vec& v, const double& vastr, double u, double& uastr
);

arma::uword gbp4d_solver_dpp_main_create_nastr(const arma::vec& p, const arma::mat& ldhw, const arma::vec& m);

arma::uword gbp4d_solver_dpp_main_create_nlmt(const arma::uword nlvl, const arma::uword nastr);

arma::vec gbp4d_solver_dpp_prep_create_p(const arma::mat& ldhw, const arma::vec& m);


// gbp4q
// @param p
//  p profit of it fit into bn <vector>
//  - cluster w, cluster max(l,d,h) and area via gbp4d_solver_dpp_main_create_p()
// @param it
//  it position and scale <matrix>
//  - x, y, z, w it position and w hold in the bin <numeric>
//  - l, d, h, w it scale along x, y, z and also w <numeric>
// @param bn
//  bn scale <matrix>
//  - l, d, h, w bn scale along x, y, z and also w <numeric>
//  - l, d, h, w in row and each col is a single bn
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
class gbp4q {

public:

  arma::vec p; arma::mat it; arma::mat bn; arma::uvec k; arma::uvec f; double o; bool ok;

  // .constructor
  gbp4q(
    arma::vec p, arma::mat it, arma::mat bn, arma::uvec k, arma::uvec f, double o, bool ok
  ) : p(p), it(it), bn(bn), k(k), f(f), o(o), ok(ok) {
  }

};


// solve gbp4d w.r.t select most preferable often smallest bin from bn list
gbp4q gbp4d_solver_dpp_filt(const arma::mat& ldhw, const arma::mat& m);

// bool assert_matrix_last_col_dominant(const arma::mat& m);

void gbp4d_solver_dpp_filt_fast(
    const arma::mat& sldhw, const arma::vec& vldhw, const arma::mat& sm, const arma::mat& vm, arma::uvec& flmt
);

void gbp4d_solver_dpp_filt_slow(
    const double vldh, const arma::rowvec& vm, arma::uvec& flmt
);

void gbp4d_solver_dpp_filt_dcol(
    const arma::mat &sm, const arma::uword icol, arma::uvec& flmt
);


#endif // __BINPACK_BINPACK4D__
