#include "gbp_u.h"
#include "gbp3d_xp.h"
// [[Rcpp::plugins(cpp11)]]

// const double tol = 0.00000001;

// //' gbp3d_xp_create_xp
// //' @description
// //'  calculate extreme point in a single bin
// //' @details
// //'  core function in gbp3d_solver_dpp
// //' @param bn
// //'  bn scale <vector>
// //'  - l, d, h bn scale along x, y, z <numeric>
// //' @param it
// //'  it position and scale <matrix>
// //'  - x, y, z it position in the bin <numeric>
// //'  - l, d, h it scale along x, y, z <numeric>
// //' @family gbp3d_xp
// //' @return xp
// //'  an updated xp
// //' @note
// //'  should make sure it kt can be fit in bin outside
// //' @note
// //'  should call gbp3d_xp_update_xp whenever possible while
// //'  direct calculate xp only after it being pushing around
// //' @export
// // [[Rcpp::export]]
arma::mat gbp3d_xp_create_xp(
  const arma::vec& bn, const arma::mat &it
) {

  // init it
  // sort it by z, y, x position - mimic fit sequence
  arma::uvec ulmt = arma::linspace<arma::uvec>(2, 0, 3);

  arma::uvec idit = sort_index_via_rows(it, ulmt);

  // init xp
  arma::mat xp = arma::zeros<arma::mat>(6, 1);

  xp(3, 0) = bn(0); xp(4, 0) = bn(1); xp(5, 0) = bn(2);

  if (it.n_cols == 0) return xp;

  // prog xp
  // fit it into bin one by one - vlmt as #.it in bin
  arma::uvec vlmt = arma::zeros<arma::uvec>(0);
  for (arma::uword i = 0; i < it.n_cols; i++) {
    if (i > 0) {
      vlmt = arma::linspace<arma::uvec>(0, i - 1, i);
    }
    gbp3d_xp_update_xp(bn, it.cols(idit(vlmt)), it.col(idit(i)), xp);
  }

  return xp;
}

// //' gbp3d_xp_update_xp
// //' @description
// //'  update extreme point in a single bin
// //' @details
// //'  core function in gbp3d_solver_dpp
// //' @param bn
// //'  bn scale <vector>
// //'  - l, d, h bn scale along x, y, z <numeric>
// //' @param it
// //'  it position and scale <matrix>
// //'  - x, y, z it position in the bin <numeric>
// //'  - l, d, h it scale along x, y, z <numeric>
// //' @param kt
// //'  kt next it placing into bin <vector>
// //'  - x, y, z kt position in the bin <numeric>
// //'  - l, d, h kt scale along x, y, z <numeric>
// //' @param xp
// //'  xp extreme point position and residual space scale <matrix>
// //'  - x, y, z xp position in the bin <numeric>
// //'  - l, d, h xp residual space scale along x, y, z <numeric>
// //' @family gbp3d_xp
// //' @return xp
// //'  an updated xp
// //' @note
// //'  should make sure it kt can be fit in bin outside
// //' @export
// // [[Rcpp::export]]
void gbp3d_xp_update_xp(
  const arma::vec& bn, const arma::mat& it, const arma::vec& kt, arma::mat& xp
) {

  // init
  if (it.n_cols == 0 && kt.size() == 0) {
    xp = arma::zeros<arma::mat>(6, 1); xp(3, 0) = bn(0); xp(4, 0) = bn(1); xp(5, 0) = bn(2); return;
  }

  // construct xp input
  // remove extreme points that is taken by or fallen into kt
  // and also update it extreme point residual space w.r.t kt
  gbp3d_xp_update_xp_ikt(it, kt, xp);

  // calculate xp output
  // calculate new extreme points from 6 projections of 3 points
  arma::mat xpUpdate(6, 6); xpUpdate.fill(arma::datum::nan);

  // maxBound: track 6 projection xp position location
  // from min x-left y-behind z-bottom to max x-right y-front z-top
  // arma::vec maxBound = arma::zeros<arma::vec>(6) - 1; // init maxBound without itBnd save computing cost
  arma::vec maxBound = arma::zeros<arma::vec>(6);

  // minBound: track 6 projection xp residual space as
  // from max x-right y-front z-top to min x-left y-behind z-bottom
  arma::mat minBound = arma::zeros<arma::mat>(3, 6);

  for (arma::uword i = 0; i < 6; i++) {
    minBound(0, i) = bn(0); // init x-right of all 6 projected extreme point
    minBound(1, i) = bn(1); // init y-front of all 6 projected extreme point
    minBound(2, i) = bn(2); // init z-top   of all 6 projected extreme point
  }

  // calculate xpUpdate x, y, z extreme point position
  gbp3d_xp_update_xp_spg(it, kt, maxBound, xpUpdate);

  // calculate xpUpdate l, d, h residual space along x, y, z
  gbp3d_xp_update_rs_spg(it, kt, minBound, xpUpdate);

  // prog xpUpdate remove nan in x, y, z, l, d, h
  arma::uvec g = arma::zeros<arma::uvec>(6);

  for (arma::uword i = 0; i < 6; i++) {
    if ((xpUpdate.col(i)).has_nan()) { g(i) = 1; }
  }

  xpUpdate = xpUpdate.cols(arma::find(g == 0));

  // join xpUpdate into xp
  xp = unique_cols(arma::join_rows(xp, xpUpdate));

  // pure xp list via remove xp with residual space == 0
  // and via remove xp dominated by other xp in the list
  gbp3d_xp_purify_xp(xp);

  // sort xp via non-decreasing order of z, y, x
  arma::uvec ulmt = arma::linspace<arma::uvec>(2, 0, 3);

  xp = sort_via_rows(xp, ulmt);

}

// //' gbp3d_xp_purify_xp
// //' @description
// //'  subroutine of gbp3d_xp_update_xp
// //' @details
// //'  cleanse xp list via remove xp with residual space == 0 and remove xp dominated by other xp in list
// //' @family gbp3d_xp
// //' @export
// // [[Rcpp::export]]
void gbp3d_xp_purify_xp(
  arma::mat& xp
) {

  // remove xp with residual space == 0
  arma::uvec g0 = arma::zeros<arma::uvec>(xp.n_cols);
  for (arma::uword i = 0; i < xp.n_cols; i++) {
    if (xp(3, i) == 0 || xp(4, i) == 0 || xp(5, i) == 0) { g0(i) = 1; }
  }
  xp = xp.cols(find(g0 == 0));

  // remove xp dominated by other xp in list
  // if x, y, z < x', y', z' and l, d, h > l', d', h' then
  // (x, y, z, l, d, h) dominant (x', y', z', l', d', h')
  // so remove (x', y', z', l', d', h') from xp list
  arma::uvec g1 = arma::zeros<arma::uvec>(xp.n_cols);
  for (arma::uword i = 0; i < xp.n_cols; i++) {
    for (arma::uword j = 0; j < xp.n_cols; j++) {
      if (i != j &&
          xp(0, i) <= xp(0, j) &&
          xp(1, i) <= xp(1, j) &&
          xp(2, i) <= xp(2, j) &&
          xp(3, i) >= xp(3, j) &&
          xp(4, i) >= xp(4, j) &&
          xp(5, i) >= xp(5, j)
      ) {
        g1(j) = 1;
      }
    }
  }
  xp = xp.cols(find(g1 == 0));

}

// //' gbp3d_xp_update_xp_ikt
// //' @description
// //'  subroutine of gbp3d_xp_update_xp
// //' @details
// //'  update current extreme point xp list based on it w.r.t fit kt into bn:
// //'  - remove extreme points that taken by kt position and remove extreme points that fallen into kt
// //'  - calculate each single extreme point residual space after fit kt into bin - relative to new kt
// //' @note
// //'  gbp3d_xp_update_xp_ikt is different to gbp3d_xp_update_rs_spg which update
// //'   extreme points introduced via kt and their residual space relative to all old it
// //' @family gbp3d_xp
// //' @export
// // [[Rcpp::export]]
void gbp3d_xp_update_xp_ikt(
  const arma::mat& it, const arma::vec& kt, arma::mat& xp
) {

  // remove extreme points that is taken by or fallen into kt
  arma::uvec vlmt = arma::zeros<arma::uvec>(xp.n_cols);

  for(arma::uword i = 0; i < xp.n_cols; i++) {
    if (kt(0) <= xp(0, i) && xp(0, i) < kt(0) + kt(3) &&
        kt(1) <= xp(1, i) && xp(1, i) < kt(1) + kt(4) &&
        kt(2) <= xp(2, i) && xp(2, i) < kt(2) + kt(5) )
      vlmt(i) = 1;
  }

  xp = xp.cols(find(vlmt == 0));

  // and also update it extreme point residual space w.r.t kt
  for (arma::uword i = 0; i < xp.n_cols; i++) {

    if (xp(0, i) <= kt(0)         &&
        xp(1, i) >= kt(1)         &&
        xp(1, i) <  kt(1) + kt(4) &&
        xp(2, i) >= kt(2)         &&
        xp(2, i) <  kt(2) + kt(5)
    ) {
      xp(3, i) = std::min(xp(3, i), kt(0) - xp(0, i));
    }

    if (xp(1, i) <= kt(1)         &&
        xp(2, i) >= kt(2)         &&
        xp(2, i) <  kt(2) + kt(5) &&
        xp(0, i) >= kt(0)         &&
        xp(0, i) <  kt(0) + kt(3)
    ) {
      xp(4, i) = std::min(xp(4, i), kt(1) - xp(1, i));
    }

    if (xp(2, i) <= kt(2)         &&
        xp(0, i) >= kt(0)         &&
        xp(0, i) <  kt(0) + kt(3) &&
        xp(1, i) >= kt(1)         &&
        xp(1, i) <  kt(1) + kt(4)
    ) {
      xp(5, i) = std::min(xp(5, i), kt(2) - xp(2, i));
    }

  }

}

// //' gbp3d_xp_update_rs_spg
// //' @description
// //'  subroutine of gbp3d_xp_update_xp
// //' @details
// //'  calculate residual space of projected kt xp over each single it in bin
// //' @family gbp3d_xp
// //' @export
// // [[Rcpp::export]]
void gbp3d_xp_update_rs_spg(
  const arma::mat& it, const arma::vec& kt,
  arma::mat& minBound, arma::mat& xpUpdate
) {

  for (arma::uword i = 0; i < it.n_cols; i++) {
    gbp3d_xp_update_minbnd(it.col(i), kt, minBound, xpUpdate);
  }

  for (arma::uword i = 0; i < 6; i++) {
    xpUpdate(3, i) = minBound(0, i) - xpUpdate(0, i);
    xpUpdate(4, i) = minBound(1, i) - xpUpdate(1, i);
    xpUpdate(5, i) = minBound(2, i) - xpUpdate(2, i);
  }

}

// //' gbp3d_xp_update_minbnd
// //' @description
// //'  subroutine of gbp3d_xp_update_xp
// //' @details
// //'  calculate residual space of projected kt xp over each single it in bin
// //' @family gbp3d_xp
// //' @export
// // [[Rcpp::export]]
void gbp3d_xp_update_minbnd(
  const arma::vec& it, const arma::vec& kt,
  arma::mat& minBound, arma::mat& xpUpdate
) {

  // projecting kt xp -> it on reverse direction w.r.t gbp3d_xp_it_pjt_kt.png
  // arma::uvec ik = gbp3d_xp_it_qjt_kt(it, kt);

  // construct virtual kt as xpUpdate with l = 0, d = 0, h = 0
  // since residual space of extreme point is related to the x, y, z of exteme point itself
  arma::vec akt(6); arma::uvec aik(6);

  for (arma::uword i = 0; i < 6; i++) {

    // init
    // init residual space without creating another itBnd and save computation cost (skip)
    // xpUpdate(3, i) = minBound(0, i) - xpUpdate(0, i);
    // xpUpdate(4, i) = minBound(1, i) - xpUpdate(1, i);
    // xpUpdate(5, i) = minBound(2, i) - xpUpdate(2, i);

    // create a virtual kt as a single point with l = 0, d = 0, h = 0
    akt(0) = xpUpdate(0, i);
    akt(1) = xpUpdate(1, i);
    akt(2) = xpUpdate(2, i);
    akt(3) = 0.00;
    akt(4) = 0.00;
    akt(5) = 0.00;

    // projecting kt xp -> it on reverse direction w.r.t gbp3d_xp_it_pjt_kt.png
    aik = gbp3d_xp_it_qjt_kt(it, akt);

    // since akt(3) = 0.00; akt(4) = 0.00; akt(5) = 0.00; =>
    // aik(0) == aik(5); aik(1) == aik(2); aik(3) == aik(4);

    // it block on the way from extreme point to x-right
    if (aik(3) && aik(4)) {
      minBound(0, i) = std::min(it(0), minBound(0, i));
    }

    if (aik(5) && aik(0)) {
      minBound(1, i) = std::min(it(1), minBound(1, i));
    }

    if (aik(1) && aik(2)) {
      minBound(2, i) = std::min(it(2), minBound(2, i));
    }

  }

}

// //' gbp3d_xp_update_xp_spg
// //' @description
// //'  subroutine of gbp3d_xp_update_xp
// //' @details
// //'  calculate extreme point position via projecting kt xp to each single it in bin
// //' @family gbp3d_xp
// //' @export
// // [[Rcpp::export]]
void gbp3d_xp_update_xp_spg(
  const arma::mat& it, const arma::vec& kt,
  arma::vec& maxBound, arma::mat& xpUpdate
) {

  for (arma::uword i = 0; i < it.n_cols; i++) {
    gbp3d_xp_update_maxbnd(it.col(i), kt, maxBound, xpUpdate);
  }

  xpUpdate(0, 0) = kt(0) + kt(3);
  xpUpdate(1, 0) = maxBound(0);
  xpUpdate(2, 0) = kt(2);

  xpUpdate(0, 1) = kt(0) + kt(3);
  xpUpdate(1, 1) = kt(1);
  xpUpdate(2, 1) = maxBound(1);

  xpUpdate(0, 2) = kt(0);
  xpUpdate(1, 2) = kt(1) + kt(4);
  xpUpdate(2, 2) = maxBound(2);

  xpUpdate(0, 3) = maxBound(3);
  xpUpdate(1, 3) = kt(1) + kt(4);
  xpUpdate(2, 3) = kt(2);

  xpUpdate(0, 4) = maxBound(4);
  xpUpdate(1, 4) = kt(1);
  xpUpdate(2, 4) = kt(2) + kt(5);

  xpUpdate(0, 5) = kt(0);
  xpUpdate(1, 5) = maxBound(5);
  xpUpdate(2, 5) = kt(2) + kt(5);

}

// //' gbp3d_xp_update_maxbnd
// //' @description
// //'  subroutine of gbp3d_xp_update_xp
// //' @details
// //'  calculate extreme point position via projecting kt xp to each single it in bin
// //' @family gbp3d_xp
// //' @export
// // [[Rcpp::export]]
void gbp3d_xp_update_maxbnd(
  const arma::vec& it, const arma::vec& kt,
  arma::vec& maxBound, arma::mat& xpUpdate
) {

  // projecting kt xp -> it along with direction w.r.t gbp3d_xp_it_pjt_kt.png
  arma::uvec ik = gbp3d_xp_it_pjt_kt(it, kt);

  // direction x-Y: kt-x-corner-move project-toward->Y
  if (ik(0) && (it(1) + it(4) > maxBound(0))) {
    maxBound(0) = it(1) + it(4);
  }

  // direction x-Z: kt-x-corner-move project-toward->Z
  if (ik(1) && (it(2) + it(5) > maxBound(1))) {
    maxBound(1) = it(2) + it(5);

  }

  // direction y-Z: kt-y-corner-move project-toward->Z
  if (ik(2) && (it(2) + it(5) > maxBound(2))) {
    maxBound(2) = it(2) + it(5);
  }

  // direction y-X: kt-y-corner-move project-toward-X
  if (ik(3) && (it(0) + it(3) > maxBound(3))) {
    maxBound(3) = it(0) + it(3);
  }

  // direction z-X: kt-z-corner-move project-toward->X
  if (ik(4) && (it(0) + it(3) > maxBound(4))) {
    maxBound(4) = it(0) + it(3);
  }

  // direction z-Y: kt-z-corner-move project-toward->Y
  if (ik(5) && (it(1) + it(4) > maxBound(5))) {
    maxBound(5) = it(1) + it(4);
  }

}

// //' gbp3d_xp_it_qjt_kt
// //' @description
// //'  can item it take projection of new extreme point xp created by next item kt on reverse direction w.r.t gbp3d_xp_it_pjt_kt.png
// //' @inheritParams gbp3d_xp_it_pjt_kt
// //' @return ik <vector>
// //'  - xY, xZ, yZ, yX, zX, zY <boolean>
// //' @family gbp3d_xp
// //' @note
// //'  xY means reverse direction of kt xp along x-axis projecting back along y-axis - see inst/img/gbp3d_xp_it_pjt_kt.png
// //' @export
// // [[Rcpp::export]]
arma::uvec gbp3d_xp_it_qjt_kt(
  const arma::vec& it, const arma::vec& kt
) {

  arma::uvec ik(6);

  // direction x-Y
  ik(0) = (
       kt(1) + kt(4) <= it(1)
    && it(0)         <= kt(0) + kt(3)
    && kt(0) + kt(3) <  it(0) + it(3)
    && it(2)         <= kt(2)
    && kt(2)         <  it(2) + it(5)
  );

  // direction x-Z
  ik(1) = (
       kt(2) + kt(5) <= it(2)
    && it(0)         <= kt(0) + kt(3)
    && kt(0) + kt(3) <  it(0) + it(3)
    && it(1)         <= kt(1)
    && kt(1)         <  it(1) + it(4)
  );

  // direction y-Z
  ik(2) = (
       kt(2) + kt(5) <= it(2)
    && it(1)         <= kt(1) + kt(4)
    && kt(1) + kt(4) <  it(1) + it(4)
    && it(0)         <= kt(0)
    && kt(0)         <  it(0) + it(3)
  );

  // direction y-X
  ik(3) = (
       kt(0) + kt(3) <= it(0)
    && it(1)         <= kt(1) + kt(4)
    && kt(1) + kt(4) <  it(1) + it(4)
    && it(2)         <= kt(2)
    && kt(2)         <  it(2) + it(5)
  );

  // direction z-X
  ik(4) = (
       kt(0) + kt(3) <= it(0)
    && it(2)         <= kt(2) + kt(5)
    && kt(2) + kt(5) <  it(2) + it(5)
    && it(1)         <= kt(1)
    && kt(1)         <  it(1) + it(4)
  );

  // direction z-Y
  ik(5) = (
       kt(1) + kt(4) <= it(1)
    && it(2)         <= kt(2) + kt(5)
    && kt(2) + kt(5) <  it(2) + it(5)
    && it(0)         <= kt(0)
    && kt(0)         <  it(0) + it(3)
  );

  return ik;
}

// //' gbp3d_xp_it_pjt_kt
// //' @description
// //'  can item it take projection of new extreme point xp created by next item kt on the one direction w.r.t gbp3d_xp_it_pjt_kt.png
// //' @param it <vector>
// //'  - x, y, z, l, d, h <numeric>
// //' @param kt <vector>
// //'  - x, y, z, l, d, h <numeric>
// //' @return ik <vector>
// //'  - xY, xZ, yZ, yX, zX, zY <boolean>
// //' @family gbp3d_xp
// //' @note
// //'  xY means the one direction of kt xp along x-axis projecting back along y-axis - see inst/img/gbp3d_xp_it_pjt_kt.png
// //' @export
// // [[Rcpp::export]]
arma::uvec gbp3d_xp_it_pjt_kt(
  const arma::vec& it, const arma::vec& kt
) {

  arma::uvec ik(6);

  // direction x-Y
  // kt raw xp along x-axis is (kt(0) + kt(3), kt(1), kt(2)) projecting back along y-axis:
  // check it Ymax is left on kt(1) and it XZ surface include point (kt(0) + kt(3), kt(2))
  ik(0) = (
       it(1) + it(4) <= kt(1)
    && it(0)         <= kt(0) + kt(3)
    && kt(0) + kt(3) <  it(0) + it(3)
    && it(2)         <= kt(2)
    && kt(2)         <  it(2) + it(5)
  );

  // direction x-Z
  ik(1) = (
       it(2) + it(5) <= kt(2)
    && it(0)         <= kt(0) + kt(3)
    && kt(0) + kt(3) <  it(0) + it(3)
    && it(1)         <= kt(1)
    && kt(1)         <  it(1) + it(4)
  );

  // direction y-Z
  ik(2) = (
       it(2) + it(5) <= kt(2)
    && it(1)         <= kt(1) + kt(4)
    && kt(1) + kt(4) <  it(1) + it(4)
    && it(0)         <= kt(0)
    && kt(0)         <  it(0) + it(3)
  );

  // direction y-X
  ik(3) = (
       it(0) + it(3) <= kt(0)
    && it(1)         <= kt(1) + kt(4)
    && kt(1) + kt(4) <  it(1) + it(4)
    && it(2)         <= kt(2)
    && kt(2)         <  it(2) + it(5)
  );

  // direction z-X
  ik(4) = (
       it(0) + it(3) <= kt(0)
    && it(2)         <= kt(2) + kt(5)
    && kt(2) + kt(5) <  it(2) + it(5)
    && it(1)         <= kt(1)
    && kt(1)         <  it(1) + it(4)
  );

  // direction z-Y
  ik(5) = (
       it(1) + it(4) <= kt(1)
    && it(2)         <= kt(2) + kt(5)
    && kt(2) + kt(5) <  it(2) + it(5)
    && it(0)         <= kt(0)
    && kt(0)         <  it(0) + it(3)
  );

  return ik;
}

