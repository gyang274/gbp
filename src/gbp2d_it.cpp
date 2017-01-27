#include "gbp2d_xp.h"
#include "gbp2d_it.h"

// gbp2d_it_create_ktlist
// @description
//  create ktlist from itlist
// @details
//  core function in gbp2d_solver_dpp
//   select highest profitable it not yet fit into bn and return all possbile fit w.r.t xp and orientation
// @param bn
//  bn scale <vector>
//  - l, d bn scale along x and y <numeric>
// @param it
//  it position and scale <matrix>
//  - x, y it position in the bin <numeric>
//  - l, d it scale along x and y <numeric>
// @param xp
//  xp extreme point position and residual space scale <matrix>
//  - x, y xp position in the bin <numeric>
//  - l, d xp residual space scale along x and y <numeric>
// @param kinit
//  kt candidate scale without position <matrix>
//  - l, d kt scale along x and y which open to orientation <numeric>
// @param nlmt
//  nlmt: limit on ktlist n max-value
// @family gbp2d_it
// @return Ktlist2d
// @note
//  should make sure it kt can be fit in bin outside
Ktlist2d gbp2d_it_create_ktlist(
  const arma::vec& bn, const arma::mat& it, const arma::mat& xp, const arma::vec& ktinit, const arma::uword nlmt
) {

  // init kt ld
  arma::mat ktldht = gbp2d_it_create_ktldht(ktinit);

  // init kt
  arma::mat kt(4, ktldht.n_cols * xp.n_cols, arma::fill::zeros);

  // ktinit -> kt w.r.t all feasible fit w.r.t orientation and xp position - impose xp residual space as fast filt
  arma::uvec vlmt = arma::zeros<arma::uvec>(ktldht.n_cols * xp.n_cols);

  arma::uword ij = 0;

  // fast filt w.r.t extreme point residual space
  for (arma::uword i = 0; i < ktldht.n_cols; i++) {
    for (arma::uword j = 0; j < xp.n_cols; j++) {
      if (ktldht(0, i) <= xp(2, j) &&
          ktldht(1, i) <= xp(3, j)
      ) {
        ij = i * xp.n_cols + j;
        kt(0, ij) = xp(0, j);
        kt(1, ij) = xp(1, j);
        kt(2, ij) = ktldht(0, i);
        kt(3, ij) = ktldht(1, i);
        vlmt(ij) = 1;
      }
    }
  }
  kt = kt.cols(arma::find(vlmt == 1));

  // ktinit -> kt w.r.t all feasible fit w.r.t orientation and xp position - impost no conflict with existing it as thorough slow filt
  vlmt = arma::ones<arma::uvec>(kt.n_cols); // re-use vlmt as no conflict indicator vector

  // slow filt w.r.t no conflict all existing it
  for (arma::uword i = 0; i < kt.n_cols; i++) {
    for (arma::uword j = 0; j < it.n_cols; j++) {
      if (!(
          kt(0, i) + kt(2, i) <= it(0, j) || it(0, j) + it(2, j) <= kt(0, i) ||
          kt(1, i) + kt(3, i) <= it(1, j) || it(1, j) + it(3, j) <= kt(1, i)
      )) {
        vlmt(i) = 0; break;
      }
    }
  }
  kt = kt.cols(arma::find(vlmt == 1));

  // init ktlist
  arma::vec s(kt.n_cols); arma::field<arma::mat> xplist(kt.n_cols);

  arma::mat xpWithKt;

  for (arma::uword i = 0; i < kt.n_cols; i++) {

    xpWithKt = xp; gbp2d_xp_update_xp(bn, it, kt.col(i), xpWithKt);

    xplist(i) = xpWithKt; s(i) = gbp2d_it_scorer_ktlist(xpWithKt);

  }

  gbp2d_it_purify_ktlist(kt, xplist, s, nlmt);

  Ktlist2d ktlist(kt.n_cols, kt, xplist, s);

  return ktlist;
}

// //' gbp2d_it_purify_ktlist
// //' @description
// //'  subroutine of gbp2d_it_create_ktlist
// //' @details
// //'  purify ktlist with s and nlmt
// //' @family gbp2d_it
// //' @export
// // [[Rcpp::export]]
void gbp2d_it_purify_ktlist(arma::mat& kt, arma::field<arma::mat>& xplist, arma::vec& s, const arma::uword nlmt) {

  // init
  arma::uvec sidx = arma::sort_index(s);

  if ((nlmt == 0) || (nlmt >= kt.n_cols)) {

    arma::field<arma::mat> xplist0(kt.n_cols);

    for (arma::uword i = 0; i < s.size(); i++) {
      xplist0(i) = xplist(sidx(i));
    }

    kt = kt.cols(sidx);  xplist = xplist0; s = s(sidx);

  } else {

    arma::uvec slmt = arma::linspace<arma::uvec>(0, nlmt - 1, nlmt);

    arma::field<arma::mat> xplist0(nlmt);

    for (arma::uword i = 0; i < nlmt; i++) {
      xplist0(i) = xplist(sidx(i));
    }

    kt = kt.cols(sidx(slmt)); xplist = xplist0; s = s(sidx(slmt));

  }

}

// //' gbp2d_it_scorer_ktlist
// //' @description
// //'  subroutine of gbp2d_it_create_ktlist
// //' @details
// //'  score kt with number of new extreme point and available scales
// //' @family gbp2d_it
// //' @export
// // [[Rcpp::export]]
double gbp2d_it_scorer_ktlist(const arma::mat& xpWithKt) {

  // init
  double s; arma::rowvec rs; // arma::uvec ulmt = arma::linspace<arma::uvec>(2, 3, 2);

  // calculcate available space via a single extreme point
  rs = xpWithKt.row(2) % xpWithKt.row(3);

  rs = rs / arma::sum(rs); // % and / element-wise mutipilication and division

  // score a xpWithKt as the entropy of the available space
  // score s is the smaller the better: prefer the configuration where less and large available dominate
  s = - arma::sum(rs % arma::log(rs));

  return s;
}

// //' gbp2d_it_create_ktldht
// //' @description
// //'  subroutine of gbp2d_it_create_ktlist
// //' @details
// //'  create all possible none duplicated orientations via mapping l d:
// //'  ktinit [x y l d] -> [l d; d l;]
// //'  also note ktinit x y are placeholder simply 0 with no meaning yet
// //' @return ktldht
// //'  all possible none duplicated orientations of ktinit in 2 x N matrix, N = 1, 2.
// //' @family gbp2d_it
// //' @export
// // [[Rcpp::export]]
arma::mat gbp2d_it_create_ktldht(const arma::vec& ktinit) {

  // init
  double l = ktinit(2);
  double d = ktinit(3);

  // main
  arma::mat ktldht(2, 2); arma::uvec vlmt = arma::zeros<arma::uvec>(2);

  // create all possible none duplicated orientations via mapping:
  // ktinit [l d] -> [l d; d l;]

  ktldht(0, 0) = l; ktldht(1, 0) = d;
  ktldht(0, 1) = d; ktldht(1, 1) = l;

  vlmt(0) = 1;

  if (l == d) {

    // vlmt(0) = 1;

  } else {

    vlmt.fill(1);

  }

  return ktldht.cols(arma::find(vlmt == 1));
}
