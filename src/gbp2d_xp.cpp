#include "gbp_u.h"
#include "gbp2d_xp.h"
// [[Rcpp::plugins(cpp11)]]

// const double tol = 0.00000001;

// //' gbp2d_xp_create_xp
// //' @description
// //'  calculate extreme point in a single bin
// //' @details
// //'  core function in gbp2d_solver_dpp
// //' @param bn
// //'  bn scale <vector>
// //'  - l, d bn scale along x and y <numeric>
// //' @param it
// //'  it position and scale <matrix>
// //'  - x, y it position in the bin <numeric>
// //'  - l, d it scale along x and y <numeric>
// //' @family gbp2d_xp
// //' @return xp
// //'  an updated xp
// //' @note
// //'  should make sure it kt can be fit in bin outside
// //' @note
// //'  should call gbp2d_xp_update_xp whenever possible while
// //'  direct calculate xp only after it being pushing around
// //' @export
// // [[Rcpp::export]]
arma::mat gbp2d_xp_create_xp(
  const arma::vec& bn, const arma::mat &it
) {

  // init it
  // sort it by y, x position - mimic fit sequence
  arma::uvec ulmt = arma::linspace<arma::uvec>(1, 0, 2);

  arma::uvec idit = sort_index_via_rows(it, ulmt);

  // init xp
  arma::mat xp = arma::zeros<arma::mat>(4, 1);

  xp(2, 0) = bn(0); xp(3, 0) = bn(1);

  if (it.n_cols == 0) return xp;

  // prog xp
  // fit it into bin one by one - vlmt as #.it in bin
  arma::uvec vlmt = arma::zeros<arma::uvec>(0);
  for (arma::uword i = 0; i < it.n_cols; i++) {
    if (i > 0) {
      vlmt = arma::linspace<arma::uvec>(0, i - 1, i);
    }
    gbp2d_xp_update_xp(bn, it.cols(idit(vlmt)), it.col(idit(i)), xp);
  }

  return xp;
}

// //' gbp2d_xp_update_xp
// //' @description
// //'  update extreme point in a single bin
// //' @details
// //'  core function in gbp2d_solver_dpp
// //' @param bn
// //'  bn scale <vector>
// //'  - l, d bn scale along x and y <numeric>
// //' @param it
// //'  it position and scale <matrix>
// //'  - x, y it position in the bin <numeric>
// //'  - l, d it scale along x and y <numeric>
// //' @param kt
// //'  kt next it placing into bin <vector>
// //'  - x, y kt position in the bin <numeric>
// //'  - l, d kt scale along x and y <numeric>
// //' @param xp
// //'  xp extreme point position and residual space scale <matrix>
// //'  - x, y xp position in the bin <numeric>
// //'  - l, d xp residual space scale along x and y <numeric>
// //' @family gbp2d_xp
// //' @return xp
// //'  an updated xp
// //' @note
// //'  should make sure it kt can be fit in bin outside
// //' @export
// // [[Rcpp::export]]
void gbp2d_xp_update_xp(
  const arma::vec& bn, const arma::mat& it, const arma::vec& kt, arma::mat& xp
) {

  // init
  if (it.n_cols == 0 && kt.size() == 0) {
    xp = arma::zeros<arma::mat>(4, 1); xp(2, 0) = bn(0); xp(3, 0) = bn(1); return;
  }

  // construct xp input
  // remove extreme points that is taken by or fallen into kt
  // and also update it extreme point residual space w.r.t kt
  gbp2d_xp_update_xp_ikt(it, kt, xp);

  // calculate xp output
  // calculate new extreme points from 2 projections of 2 points
  arma::mat xpUpdate(4, 2); xpUpdate.fill(arma::datum::nan);

  // maxBound: track 2 projection xp position location
  // from min x-left y-bottom to max x-right y-top
  arma::vec maxBound = arma::zeros<arma::vec>(2);
  // minBound: track 2 projection xp residual space as
  // from max x-right y-top to min x-left y-bottom
  arma::mat minBound = arma::zeros<arma::mat>(2, 2);
  for (arma::uword i = 0; i < 2; i++) {
    minBound(0, i) = bn(0); // init x-right of all 2 projected extreme point
    minBound(1, i) = bn(1); // init y-top   of all 2 projected extreme point
  }

  // calculate xpUpdate x, y extreme point position
  gbp2d_xp_update_xp_spg(it, kt, maxBound, xpUpdate);

  // calculate xpUpdate l, d residual space along x, y
  gbp2d_xp_update_rs_spg(it, kt, minBound, xpUpdate);

  // prog xpUpdate remove nan in x, y, l, d
  arma::uvec g = arma::zeros<arma::uvec>(2);

  for (arma::uword i = 0; i < 2; i++) {
    if ((xpUpdate.col(i)).has_nan()) { g(i) = 1; }
  }

  xpUpdate = xpUpdate.cols(arma::find(g == 0));

  // join xpUpdate into xp
  xp = unique_cols(arma::join_rows(xp, xpUpdate));

  // pure xp list via remove xp with residual space == 0
  // and via remove xp dominated by other xp in the list
  gbp2d_xp_purify_xp(xp);

  // sort xp via non-decreasing order of y, x
  arma::uvec ulmt = arma::linspace<arma::uvec>(1, 0, 2);

  xp = sort_via_rows(xp, ulmt);

}

// //' gbp2d_xp_purify_xp
// //' @description
// //'  subroutine of gbp2d_xp_update_xp
// //' @details
// //'  cleanse xp list via remove xp with residual space == 0 and remove xp dominated by other xp in list
// //' @family gbp2d_xp
// //' @export
// // [[Rcpp::export]]
void gbp2d_xp_purify_xp(
  arma::mat& xp
) {

  // remove xp with residual space == 0
  arma::uvec g0 = arma::zeros<arma::uvec>(xp.n_cols);
  for (arma::uword i = 0; i < xp.n_cols; i++) {
    if (xp(2, i) == 0 || xp(3, i) == 0) { g0(i) = 1; }
  }
  xp = xp.cols(find(g0 == 0));

  // remove xp dominated by other xp in list
  // if x, y < x', y' and l, d > l', d' then
  // (x, y, l, d) dominant (x', y', l', d')
  // so remove (x', y', l', d') from xp list
  arma::uvec g1 = arma::zeros<arma::uvec>(xp.n_cols);
  for (arma::uword i = 0; i < xp.n_cols; i++) {
    for (arma::uword j = 0; j < xp.n_cols; j++) {
      if (i != j &&
          xp(0, i) <= xp(0, j) &&
          xp(1, i) <= xp(1, j) &&
          xp(2, i) >= xp(2, j) &&
          xp(3, i) >= xp(3, j)
      ) {
        g1(j) = 1;
      }
    }
  }
  xp = xp.cols(find(g1 == 0));

}

// //' gbp2d_xp_update_xp_ikt
// //' @description
// //'  subroutine of gbp2d_xp_update_xp
// //' @details
// //'  update current extreme point xp list based on it w.r.t fit kt into bn:
// //'  - remove extreme points that taken by kt position and remove extreme points that fallen into kt
// //'  - calculate each single extreme point residual space after fit kt into bin - relative to new kt
// //' @note
// //'  gbp2d_xp_update_xp_ikt is different to gbp2d_xp_update_rs_spg which update
// //'   extreme points introduced via kt and their residual space relative to all old it
// //' @family gbp2d_xp
// //' @export
// // [[Rcpp::export]]
void gbp2d_xp_update_xp_ikt(
  const arma::mat& it, const arma::vec& kt, arma::mat& xp
) {

  // remove extreme points that is taken by or fallen into kt
  arma::uvec vlmt = arma::zeros<arma::uvec>(xp.n_cols);

  for(arma::uword i = 0; i < xp.n_cols; i++) {
    if (kt(0) <= xp(0, i) && xp(0, i) < kt(0) + kt(2) &&
        kt(1) <= xp(1, i) && xp(1, i) < kt(1) + kt(3) )
      vlmt(i) = 1;
  }

  xp = xp.cols(find(vlmt == 0));

  // and also update it extreme point residual space w.r.t kt
  for (arma::uword i = 0; i < xp.n_cols; i++) {

    if (xp(0, i) <= kt(0)         &&
        xp(1, i) >= kt(1)         &&
        xp(1, i) <  kt(1) + kt(3)
    ) {
      xp(2, i) = std::min(xp(2, i), kt(0) - xp(0, i));
    }

    if (xp(1, i) <= kt(1)         &&
        xp(0, i) >= kt(0)         &&
        xp(0, i) <  kt(0) + kt(2)
    ) {
      xp(3, i) = std::min(xp(3, i), kt(1) - xp(1, i));
    }

  }

}

// //' gbp2d_xp_update_rs_spg
// //' @description
// //'  subroutine of gbp2d_xp_update_xp
// //' @details
// //'  calculate residual space of projected kt xp over each single it in bin
// //' @family gbp2d_xp
// //' @export
// // [[Rcpp::export]]
void gbp2d_xp_update_rs_spg(
  const arma::mat& it, const arma::vec& kt,
  arma::mat& minBound, arma::mat& xpUpdate
) {

  for (arma::uword i = 0; i < it.n_cols; i++) {
    gbp2d_xp_update_minbnd(it.col(i), kt, minBound, xpUpdate);
  }

  for (arma::uword i = 0; i < 2; i++) {
    xpUpdate(2, i) = minBound(0, i) - xpUpdate(0, i);
    xpUpdate(3, i) = minBound(1, i) - xpUpdate(1, i);
  }

}

// //' gbp2d_xp_update_minbnd
// //' @description
// //'  subroutine of gbp2d_xp_update_xp
// //' @details
// //'  calculate residual space of projected kt xp over each single it in bin
// //' @family gbp2d_xp
// //' @export
// // [[Rcpp::export]]
void gbp2d_xp_update_minbnd(
  const arma::vec& it, const arma::vec& kt,
  arma::mat& minBound, arma::mat& xpUpdate
) {

  // projecting kt xp -> it on reverse direction w.r.t gbp2d_xp_it_pjt_kt.png
  // arma::uvec ik = gbp2d_xp_it_qjt_kt(it, kt);

  // construct virtual kt as xpUpdate with l = 0, d = 0
  // since residual space of extreme point is related to the x, y of exteme point itself
  arma::vec akt(4); arma::uvec aik(2);

  for (arma::uword i = 0; i < 2; i++) {

    // init
    // init residual space without creating another itBnd and save computation cost (skip)
    // xpUpdate(2, i) = minBound(0, i) - xpUpdate(0, i);
    // xpUpdate(3, i) = minBound(1, i) - xpUpdate(1, i);

    // create a virtual kt as a single point with l = 0, d = 0
    akt(0) = xpUpdate(0, i);
    akt(1) = xpUpdate(1, i);
    akt(2) = 0.00;
    akt(3) = 0.00;

    // projecting kt xp -> it on reverse direction w.r.t gbp2d_xp_it_pjt_kt.png
    aik = gbp2d_xp_it_qjt_kt(it, akt);

    // it block on the way from extreme point to x-right
    if (aik(1)) {
      minBound(0, i) = std::min(it(0), minBound(0, i));
    }

    // it block on the way from extreme point to y-upper
    if (aik(0)) {
      minBound(1, i) = std::min(it(1), minBound(1, i));
    }

  }

}

// //' gbp2d_xp_update_xp_spg
// //' @description
// //'  subroutine of gbp2d_xp_update_xp
// //' @details
// //'  calculate extreme point position via projecting kt xp to each single it in bin
// //' @family gbp2d_xp
// //' @export
// // [[Rcpp::export]]
void gbp2d_xp_update_xp_spg(
  const arma::mat& it, const arma::vec& kt,
  arma::vec& maxBound, arma::mat& xpUpdate
) {

  for (arma::uword i = 0; i < it.n_cols; i++) {
    gbp2d_xp_update_maxbnd(it.col(i), kt, maxBound, xpUpdate);
  }

  xpUpdate(0, 0) = kt(0) + kt(2);
  xpUpdate(1, 0) = maxBound(0);

  xpUpdate(0, 1) = maxBound(1);
  xpUpdate(1, 1) = kt(1) + kt(3);

}

// //' gbp2d_xp_update_maxbnd
// //' @description
// //'  subroutine of gbp2d_xp_update_xp
// //' @details
// //'  calculate extreme point position via projecting kt xp to each single it in bin
// //' @family gbp2d_xp
// //' @export
// // [[Rcpp::export]]
void gbp2d_xp_update_maxbnd(
  const arma::vec& it, const arma::vec& kt,
  arma::vec& maxBound, arma::mat& xpUpdate
) {

  // projecting kt xp -> it along with direction w.r.t gbp2d_xp_it_pjt_kt.png
  arma::uvec ik = gbp2d_xp_it_pjt_kt(it, kt);

  // direction x-Y: kt-x-corner-move project-toward->Y
  if (ik(0) && (it(1) + it(3) > maxBound(0))) {
    maxBound(0) = it(1) + it(3);
  }

  // direction y-X: kt-y-corner-move project-toward-X
  if (ik(1) && (it(0) + it(2) > maxBound(1))) {
    maxBound(1) = it(0) + it(2);
  }

}

// //' gbp2d_xp_it_qjt_kt
// //' @description
// //'  can item it take projection of new extreme point xp created by next item kt on reverse direction w.r.t gbp2d_xp_it_pjt_kt.png
// //' @inheritParams gbp2d_xp_it_pjt_kt
// //' @return ik <vector>
// //'  - xY, yX <boolean>
// //' @family gbp2d_xp
// //' @note
// //'  xY means reverse direction of kt xp along x-axis projecting back along y-axis - see inst/img/gbp2d_xp_it_pjt_kt.png
// //' @export
// // [[Rcpp::export]]
arma::uvec gbp2d_xp_it_qjt_kt(
  const arma::vec& it, const arma::vec& kt
) {

  arma::uvec ik(2);

  // direction x-Y
  ik(0) = (
       kt(1) + kt(3) <= it(1)
    && it(0)         <= kt(0) + kt(2)
    && kt(0) + kt(2) <  it(0) + it(2)
  );

  // direction y-X
  ik(1) = (
       kt(0) + kt(2) <= it(0)
    && it(1)         <= kt(1) + kt(3)
    && kt(1) + kt(3) <  it(1) + it(3)
  );

  return ik;
}

// //' gbp2d_xp_it_pjt_kt
// //' @description
// //'  can item it take projection of new extreme point xp created by next item kt on the one direction w.r.t gbp2d_xp_it_pjt_kt.png
// //' @param it <vector>
// //'  - x, y, l, d <numeric>
// //' @param kt <vector>
// //'  - x, y, l, d <numeric>
// //' @return ik <vector>
// //'  - xY, yX <boolean>
// //' @family gbp2d_xp
// //' @note
// //'  xY means the one direction of kt xp along x-axis projecting back along y-axis - see inst/img/gbp2d_xp_it_pjt_kt.png
// //' @export
// // [[Rcpp::export]]
arma::uvec gbp2d_xp_it_pjt_kt(
  const arma::vec& it, const arma::vec& kt
) {

  arma::uvec ik(2);

  // direction x-Y
  // kt raw xp along x-axis is (kt(0) + kt(2), kt(1)) projecting back along y-axis:
  // check it Ymax is under on kt(1) and it X span area include point kt(0) + kt(2)
  ik(0) = (
       it(1) + it(3) <= kt(1)
    && it(0)         <= kt(0) + kt(2)
    && kt(0) + kt(2) <  it(0) + it(2)
  );

  // direction y-X
  ik(1) = (
       it(0) + it(2) <= kt(0)
    && it(1)         <= kt(1) + kt(3)
    && kt(1) + kt(3) <  it(1) + it(3)
  );

  return ik;
}

