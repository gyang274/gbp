#include "gbp4d_xp.h"
#include "gbp4d_it.h"

// gbp4d_it_create_ktlist
// @description
//  create ktlist from itlist
// @details
//  core function in gbp4d_solver_dpp
//   select highest profitable it not yet fit into bn and return all possbile fit w.r.t xp and orientation
// @param bn
//  bn scale <vector>
//  - l, d, h, w bn scale along x, y, z and w <numeric>
// @param it
//  it position and scale <matrix>
//  - x, y, z, w it position and w in the bin <numeric>
//  - l, d, h, w it scale along x, y, z and w <numeric>
// @param xp
//  xp extreme point position and residual space scale <matrix>
//  - x, y, z, w xp position and w in the bin <numeric>
//  - l, d, h, w xp residual space scale along x, y, z and w <numeric>
// @param kinit
//  kt candidate scale without position <matrix>
//  - l, d, h, w kt scale along x, y, z, and w which open to orientation <numeric>
// @param nlmt
//  nlmt: limit on ktlist n max-value
// @family gbp4d_it
// @return Ktlist4d
// @note
//  should make sure it kt can be fit in bin outside
Ktlist4d gbp4d_it_create_ktlist(
  const arma::vec& bn, const arma::mat& it, const arma::mat& xp, const arma::vec& ktinit, const arma::uword nlmt
) {

  // init kt ldh
  arma::mat ktldht = gbp4d_it_create_ktldht(ktinit);

  // init kt
  arma::mat kt(8, ktldht.n_cols * xp.n_cols, arma::fill::zeros);

  // ktinit -> kt w.r.t all feasible fit w.r.t orientation and xp position - impose xp residual space as fast filt
  arma::uvec vlmt = arma::zeros<arma::uvec>(ktldht.n_cols * xp.n_cols);

  arma::uword ij = 0;

  // fast filt w.r.t extreme point residual space
  for (arma::uword i = 0; i < ktldht.n_cols; i++) {
    for (arma::uword j = 0; j < xp.n_cols; j++) {
      if (ktldht(0, i) <= xp(4, j) &&
          ktldht(1, i) <= xp(5, j) &&
          ktldht(2, i) <= xp(6, j) &&
          ktldht(3, i) <= xp(7, j)
      ) {
        ij = i * xp.n_cols + j;
        kt(0, ij) = xp(0, j);
        kt(1, ij) = xp(1, j);
        kt(2, ij) = xp(2, j);
        kt(3, ij) = xp(3, j);
        kt(4, ij) = ktldht(0, i);
        kt(5, ij) = ktldht(1, i);
        kt(6, ij) = ktldht(2, i);
        kt(7, ij) = ktldht(3, i);
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
        kt(0, i) + kt(4, i) <= it(0, j) || it(0, j) + it(4, j) <= kt(0, i) ||
        kt(1, i) + kt(5, i) <= it(1, j) || it(1, j) + it(5, j) <= kt(1, i) ||
        kt(2, i) + kt(6, i) <= it(2, j) || it(2, j) + it(6, j) <= kt(2, i)
      // kt(3, i) == sum(it(7, j)) // weight on separate single dimension
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

    xpWithKt = xp; gbp4d_xp_update_xp(bn, it, kt.col(i), xpWithKt);

    xplist(i) = xpWithKt; s(i) = gbp4d_it_scorer_ktlist(xpWithKt);

  }

  gbp4d_it_purify_ktlist(kt, xplist, s, nlmt);

  Ktlist4d ktlist(kt.n_cols, kt, xplist, s);

  return ktlist;
}

// //' gbp4d_it_purify_ktlist
// //' @description
// //'  subroutine of gbp4d_it_create_ktlist
// //' @details
// //'  purify ktlist with s and nlmt
// //' @family gbp4d_it
// //' @export
// // [[Rcpp::export]]
void gbp4d_it_purify_ktlist(arma::mat& kt, arma::field<arma::mat>& xplist, arma::vec& s, const arma::uword nlmt) {//

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

// //' gbp4d_it_scorer_ktlist
// //' @description
// //'  subroutine of gbp4d_it_create_ktlist
// //' @details
// //'  score kt with number of new extreme point and available scales
// //' @family gbp4d_it
// //' @export
// // [[Rcpp::export]]
double gbp4d_it_scorer_ktlist(const arma::mat& xpWithKt) {

  // init
  double s; arma::rowvec rs; // arma::uvec ulmt = arma::linspace<arma::uvec>(4, 7, 4);

  // calculcate available space via a single extreme point
  rs = xpWithKt.row(4) % xpWithKt.row(5) % xpWithKt.row(6);

  rs = rs / arma::sum(rs); // % and / element-wise mutipilication and division

  // score a xpWithKt as the entropy of the available space
  // score s is the smaller the better: prefer the configuration where less and large available dominate
  s = - arma::sum(rs % arma::log(rs));

  return s;
}

// //' gbp4d_it_create_ktldht
// //' @description
// //'  subroutine of gbp4d_it_create_ktlist
// //' @details
// //'  create all possible none duplicated orientations via mapping l d h w:
// //'  ktinit [x y z w l d h w] -> [l d h w; l h d w; d l h w; d h l w; h l d w; h d l w;]
// //'  also note ktinit x y z w are placeholder simply 0 with no meaning yet
// //' @return ktldht
// //'  all possible none duplicated orientations of ktinit in 4 x N matrix, N = 1, 3, 6.
// //' @family gbp4d_it
// //' @export
// // [[Rcpp::export]]
arma::mat gbp4d_it_create_ktldht(const arma::vec& ktinit) {

  // init
  double l = ktinit(4);
  double d = ktinit(5);
  double h = ktinit(6);
  double w = ktinit(7);

  // main
  arma::mat ktldht(4, 6); arma::uvec vlmt = arma::zeros<arma::uvec>(6);

  // create all possible none duplicated orientations via mapping:
  // ktinit [l d h w] -> [l d h w; l h d w; d l h w; d h l w; h l d w; h d l w;]

  ktldht(0, 0) = l; ktldht(1, 0) = d; ktldht(2, 0) = h; ktldht(3, 0) = w;
  ktldht(0, 1) = l; ktldht(1, 1) = h; ktldht(2, 1) = d; ktldht(3, 1) = w;
  ktldht(0, 2) = d; ktldht(1, 2) = l; ktldht(2, 2) = h; ktldht(3, 2) = w;
  ktldht(0, 3) = d; ktldht(1, 3) = h; ktldht(2, 3) = l; ktldht(3, 3) = w;
  ktldht(0, 4) = h; ktldht(1, 4) = l; ktldht(2, 4) = d; ktldht(3, 4) = w;
  ktldht(0, 5) = h; ktldht(1, 5) = d; ktldht(2, 5) = l; ktldht(3, 5) = w;

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
