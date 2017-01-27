#include "gbp3d_xp.h"
#include "gbp3d_it.h"

// gbp3d_it_create_ktlist
// @description
//  create ktlist from itlist
// @details
//  core function in gbp3d_solver_dpp
//   select highest profitable it not yet fit into bn and return all possbile fit w.r.t xp and orientation
// @param bn
//  bn scale <vector>
//  - l, d, h bn scale along x, y, z <numeric>
// @param it
//  it position and scale <matrix>
//  - x, y, z it position in the bin <numeric>
//  - l, d, h it scale along x, y, z <numeric>
// @param xp
//  xp extreme point position and residual space scale <matrix>
//  - x, y, z xp position in the bin <numeric>
//  - l, d, h xp residual space scale along x, y, z <numeric>
// @param kinit
//  kt candidate scale without position <matrix>
//  - l, d, h kt scale along x, y, z which open to orientation <numeric>
// @param nlmt
//  nlmt: limit on ktlist n max-value
// @family gbp3d_it
// @return Ktlist3d
// @note
//  should make sure it kt can be fit in bin outside
Ktlist3d gbp3d_it_create_ktlist(
  const arma::vec& bn, const arma::mat& it, const arma::mat& xp, const arma::vec& ktinit, const arma::uword nlmt
) {

  // init kt ldh
  arma::mat ktldht = gbp3d_it_create_ktldht(ktinit);

  // init kt
  arma::mat kt(6, ktldht.n_cols * xp.n_cols, arma::fill::zeros);

  // ktinit -> kt w.r.t all feasible fit w.r.t orientation and xp position - impose xp residual space as fast filt
  arma::uvec vlmt = arma::zeros<arma::uvec>(ktldht.n_cols * xp.n_cols);

  arma::uword ij = 0;

  // fast filt w.r.t extreme point residual space
  for (arma::uword i = 0; i < ktldht.n_cols; i++) {
    for (arma::uword j = 0; j < xp.n_cols; j++) {
      if (ktldht(0, i) <= xp(3, j) &&
          ktldht(1, i) <= xp(4, j) &&
          ktldht(2, i) <= xp(5, j)
      ) {
        ij = i * xp.n_cols + j;
        kt(0, ij) = xp(0, j);
        kt(1, ij) = xp(1, j);
        kt(2, ij) = xp(2, j);
        kt(3, ij) = ktldht(0, i);
        kt(4, ij) = ktldht(1, i);
        kt(5, ij) = ktldht(2, i);
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
        kt(0, i) + kt(3, i) <= it(0, j) || it(0, j) + it(3, j) <= kt(0, i) ||
        kt(1, i) + kt(4, i) <= it(1, j) || it(1, j) + it(4, j) <= kt(1, i) ||
        kt(2, i) + kt(5, i) <= it(2, j) || it(2, j) + it(5, j) <= kt(2, i)
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

    xpWithKt = xp; gbp3d_xp_update_xp(bn, it, kt.col(i), xpWithKt);

    xplist(i) = xpWithKt; s(i) = gbp3d_it_scorer_ktlist(xpWithKt);

  }

  gbp3d_it_purify_ktlist(kt, xplist, s, nlmt);

  Ktlist3d ktlist(kt.n_cols, kt, xplist, s);

  return ktlist;
}

// //' gbp3d_it_purify_ktlist
// //' @description
// //'  subroutine of gbp3d_it_create_ktlist
// //' @details
// //'  purify ktlist with s and nlmt
// //' @family gbp3d_it
// //' @export
// // [[Rcpp::export]]
void gbp3d_it_purify_ktlist(arma::mat& kt, arma::field<arma::mat>& xplist, arma::vec& s, const arma::uword nlmt) {

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

// //' gbp3d_it_scorer_ktlist
// //' @description
// //'  subroutine of gbp3d_it_create_ktlist
// //' @details
// //'  score kt with number of new extreme point and available scales
// //' @family gbp3d_it
// //' @export
// // [[Rcpp::export]]
double gbp3d_it_scorer_ktlist(const arma::mat& xpWithKt) {

  // init
  // arma::vec s(5); arma::uvec ulmt = arma::linspace<arma::uvec>(3, 5, 3);

  // number of new xp introduced - the smaller the better
  // s(0) = xpWithKt.n_cols - xp.n_cols;

  // lost of scale on x, y, z dimensions - the smaller the better
  // s(1) = xpWithKt.row(0).max() - xp.row(0).max() + xpWithKt.row(1).max() - xp.row(1).max() + xpWithKt.row(2).max() - xp.row(2).max();

  // maximum of available space via a single extreme point - the larger the better
  // s(2) = (xpWithKt.row(3) % xpWithKt.row(4) % xpWithKt.row(5)).max(); // %: schur product: element-wise multiplication of two objects
  // maximum of available space via a single extreme point - consider negative so the smaller the better
  // s(2) = - (xpWithKt.row(3) % xpWithKt.row(4) % xpWithKt.row(5)).max(); // %: schur product: element-wise multiplication of two objects

  // maximum of available scale via a single extreme point - the larger the better
  // s(3) = xpWithKt.rows(ulmt).max();
  // maximum of available scale via a single extreme point - consider negative so the smaller the better
  // s(3) = - xpWithKt.rows(ulmt).max();

  // minimum of available scale via a single extreme point - the larger the better
  // s(4) = xpWithKt.rows(ulmt).min();
  // minimum of available scale via a single extreme point - consider negative so the smaller the better
  // s(4) = - xpWithKt.rows(ulmt).min();

  // init
  double s; arma::rowvec rs; // arma::uvec ulmt = arma::linspace<arma::uvec>(3, 5, 3);

  // calculcate available space via a single extreme point
  rs = xpWithKt.row(3) % xpWithKt.row(4) % xpWithKt.row(5);

  rs = rs / arma::sum(rs); // % and / element-wise mutipilication and division

  // score a xpWithKt as the entropy of the available space
  // score s is the smaller the better: prefer the configuration where less and large available dominate
  s = - arma::sum(rs % arma::log(rs));

  return s;
}

// //' gbp3d_it_create_ktldht
// //' @description
// //'  subroutine of gbp3d_it_create_ktlist
// //' @details
// //'  create all possible none duplicated orientations via mapping l d h:
// //'  ktinit [x y z l d h] -> [l d h; l h d; d l h; d h l; h l d; h d l;]
// //'  also note ktinit x y z are placeholder simply 0 with no meaning yet
// //' @return ktldht
// //'  all possible none duplicated orientations of ktinit in 3 x N matrix, N = 1, 3, 6.
// //' @family gbp3d_it
// //' @export
// // [[Rcpp::export]]
arma::mat gbp3d_it_create_ktldht(const arma::vec& ktinit) {

  // init
  double l = ktinit(3);
  double d = ktinit(4);
  double h = ktinit(5);

  // main
  arma::mat ktldht(3, 6); arma::uvec vlmt = arma::zeros<arma::uvec>(6);

  // create all possible none duplicated orientations via mapping:
  // ktinit [l d h] -> [l d h; l h d; d l h; d h l; h l d; h d l;]

  ktldht(0, 0) = l; ktldht(1, 0) = d; ktldht(2, 0) = h;
  ktldht(0, 1) = l; ktldht(1, 1) = h; ktldht(2, 1) = d;
  ktldht(0, 2) = d; ktldht(1, 2) = l; ktldht(2, 2) = h;
  ktldht(0, 3) = d; ktldht(1, 3) = h; ktldht(2, 3) = l;
  ktldht(0, 4) = h; ktldht(1, 4) = l; ktldht(2, 4) = d;
  ktldht(0, 5) = h; ktldht(1, 5) = d; ktldht(2, 5) = l;

  vlmt(0) = 1;

  if ((l == d) && (d == h)) {

    // vlmt(0) = 1;

  } else if (l == d) {

    // if l == d => ktldht 02 13 45 are groups within which orientation equivalent with each other

    // vlmt(0) = 1;

    vlmt(1) = 1; vlmt(4) = 1;

  } else if (d == h) {

    // if d == h => ktldht 01 24 35 are groups within which orientation equivalent with each other

    // vlmt(0) = 1;

    vlmt(2) = 1; vlmt(3) = 1;

  } else if (h == l) {

    // if h == l => ktldht 05 14 23 are groups within which orientation equivalent with each other

    // vlmt(0) = 1;

    vlmt(1) = 1; vlmt(2) = 1;

  } else {

    vlmt.fill(1);

  }

  return ktldht.cols(arma::find(vlmt == 1));
}
