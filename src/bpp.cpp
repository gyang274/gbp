#include "gbp_u.h"
#include "gbp1d.h"
#include "gbp2d.h"
#include "gbp3d.h"
#include "gbp4d.h"
#include "bpp.h"

//' bpp_solver_dpp_wrapper
//'
//' @description
//'
//'  a wrapper over bpp_solver_dpp and expose an nicer r interface
//'
//' @param it <data.frame>
//'
//'  it order itemSKU list
//'
//'  - oid: order id <integer>
//'
//'  - sku: stock keeping unit - it id <character>
//'
//'  - l, d, h, w it scale along x, y, z and w <numeric>
//'
//'  - w will be used as constraint while l, d, h will be used as both constraint and objective
//'
//'  it must be sorted w.r.t oid
//'
//' @param bn <data.frame>
//'
//'  bn a bin list
//'
//'  - id: bin id <character>
//'
//'  - l, d, h, w bn scale along x, y, z and w <numeric>
//'
//'  bn must be sorted w.r.t preference and have l >= d >= h
//'
//' @return sn <list>
//'
//'  sn solution - it order itemSKU list with tid, bid, and x, y, z <data.frame>
//'
//'  - oid: order id inherited from it <character>
//'
//'  - tid: ticket id implied one order can be packed using several ticket id <character>
//'
//'    each ticket id corresponding to a bid bin id which indicates which bin to use for packing
//'
//'  - bid: bin id which bn in bn list should be used in pakcing <character>
//'
//'  - sku: stock keeping unit it id <character>
//'
//'  - x, y, z it position in the bin <numeric>
//'
//'  - l, d, h it scale along x, y, z <numeric>
//'
//'    l, d, h is not inherited from it as it can be rotated to different orientation for packing
//'
//'  - w it weight scale inherited from it <numeric>
//'
//' @family bpp_solver_dpp
//' @export
// [[Rcpp::export]]
DataFrame bpp_solver_dpp_wrapper(
    DataFrame it, DataFrame bn
) {

  // init
  if (it.nrows() == 0) {
    return Rcpp::DataFrame::create(
      Rcpp::Named("oid") = arma::uvec(),
      Rcpp::Named("tid") = arma::uvec(),
      Rcpp::Named("bid") = arma::uvec(),
      Rcpp::Named("sku") = arma::uvec(),
      Rcpp::Named("x") = arma::vec(),
      Rcpp::Named("y") = arma::vec(),
      Rcpp::Named("z") = arma::vec(),
      Rcpp::Named("l") = arma::vec(),
      Rcpp::Named("d") = arma::vec(),
      Rcpp::Named("h") = arma::vec(),
      Rcpp::Named("w") = arma::uvec(),
      Rcpp::Named("stringsAsFactors") = false
    );
  }

  // init it

  arma::uvec id = as<arma::uvec>(it["oid"]);

  CharacterVector itsku = it["sku"];

  arma::mat ldhw(4, it.nrows(), arma::fill::zeros);

  ldhw.row(0) = as<arma::rowvec>(it["l"]);
  ldhw.row(1) = as<arma::rowvec>(it["d"]);
  ldhw.row(2) = as<arma::rowvec>(it["h"]);
  ldhw.row(3) = as<arma::rowvec>(it["w"]);

  // init bn
  CharacterVector bnidlist = bn["id"];

  arma::mat m(4, bn.nrows());

  m.row(0) = as<arma::rowvec>(bn["l"]);
  m.row(1) = as<arma::rowvec>(bn["d"]);
  m.row(2) = as<arma::rowvec>(bn["h"]);
  m.row(3) = as<arma::rowvec>(bn["w"]);

  arma::mat sm(4, bn.nrows());

  sm.rows(0, 2) = arma::sort(m.rows(0, 2), "descend", 0); sm.row(3) = m.row(3);

  m = sm; // protect: m must have l >= d >= h;

  // main
  bppSgl afit = bpp_solver_dpp(id, ldhw, m);

  arma::vec snitx = arma::trans(afit.it.row(0));
  arma::vec snity = arma::trans(afit.it.row(1));
  arma::vec snitz = arma::trans(afit.it.row(2));
  // arma::vec snitw = arma::trans(afit.it.row(3)); // w hold in bin when fit it
  arma::vec snitl = arma::trans(afit.it.row(4));
  arma::vec snitd = arma::trans(afit.it.row(5));
  arma::vec snith = arma::trans(afit.it.row(6));
  arma::vec snitw = arma::trans(afit.it.row(7)); // w of it

  CharacterVector snbid(it.nrows());

  for (int i = 0; i < it.nrows(); i++) { snbid(i) = bnidlist(afit.kb(i)); }

  return Rcpp::DataFrame::create(
    Rcpp::Named("oid") = id,
    Rcpp::Named("tid") = afit.k,
    Rcpp::Named("bid") = snbid,
    Rcpp::Named("sku") = itsku,
    Rcpp::Named("x") = snitx, // Rcpp::wrap(afit.it.row(0)), // cannot wrap subview?
    Rcpp::Named("y") = snity, // Rcpp::wrap(afit.it.row(1)), // cannot wrap subview?
    Rcpp::Named("z") = snitz, // Rcpp::wrap(afit.it.row(2)), // cannot wrap subview?
    Rcpp::Named("l") = snitl, // Rcpp::wrap(afit.it.row(4)), // cannot wrap subview?
    Rcpp::Named("d") = snitd, // Rcpp::wrap(afit.it.row(5)), // cannot wrap subview?
    Rcpp::Named("h") = snith, // Rcpp::wrap(afit.it.row(6)), // cannot wrap subview?
    Rcpp::Named("w") = snitw, // Rcpp::wrap(afit.it.row(7)), // cannot wrap subview?
    Rcpp::Named("stringsAsFactors") = false
  );

}


// bpp_solver_dpp
// @description
//  main solver of e-commerce warehouse packing algorithm
// @details
//  bpp init a list of order on sku in data.frame it - oid, sku, l, d, h, w:
//   order id oid, stock keeping unit sku, length l, depth d, height h and weight w.
//  and also a list of available bn in data.frame bn - id, l, d, h:
//   bn id, length l, depth d, height h, sorted by peference ealier smaller prefered
//  and a single weight limit wlmt applied on all bin.
//  bpp solver would solve
//   select least number of bn for packing each order w.r.t bn size and weight limit
//    and make sure the bn selected are as small as possible.
// @param id <vector>
//  id order id <integer> vector - should sorted or at least grouped w.r.t order id
// @param ldhw <matrix>
//  it order list
//  - l, d, h, w it scale along x, y, z and also w <numeric>
//  it columns should corresponding to olist
// @param m <matrix>
//  m a bin list
//  - l, d, h, w bn scale along x, y, z and also w <numeric>
//  m should sorted w.r.t preference
// @return bppSgl
// @family bpp_solver_dpp
// @export
bppSgl bpp_solver_dpp(
  const arma::uvec& id, const arma::mat& ldhw, const arma::mat& m
) {

  // init
  const arma::mat bn = m; // bn

  // if (ldhw.n_cols == 0) { return bppSgl(arma::uvec(), arma::mat(8, 0), bn, arma::uvec(), arma::uvec(), true); } // won't enter while anyway?

  // init
  arma::mat it(8, ldhw.n_cols, arma::fill::zeros); // it - x, y, z, w (w hold in bn when fit it), l, d, h, w (w of it) on row 0, 1, 2, 3, 4, 5, 6, 7

  arma::uvec xyzw_ulmt = arma::linspace<arma::uvec>(0, 3, 4);

  arma::uvec ldhw_ulmt = arma::linspace<arma::uvec>(4, 7, 4);

  it.rows(xyzw_ulmt).fill(-1); it.rows(ldhw_ulmt) = ldhw; // it - init x, y, z, w, l, d, h, w into it

  // init
  arma::uvec k(ldhw.n_cols, arma::fill::zeros); // k - ticket id vector

  // init
  arma::uvec kd(ldhw.n_cols, arma::fill::zeros); // kd - ticket bn id vector

  bool ok = true;

  // main
  arma::uword olock = id(0); arma::uword oinit = 0; arma::uword otail = 0;

  bppSgl afit(arma::uvec(), arma::mat(8, 0), arma::mat(4, 0), arma::uvec(), arma::uvec(), false);

  for (arma::uword i = 0; i <= id.size(); i++) {

    if ((i == id.size()) || (id(i) != olock)) {

      otail = i - 1;

      Rcpp::Rcout << "bpp_solver_dpp: processing order id: " << olock << " on index: " << oinit << " - " << otail << " .. "<< std::endl;

      afit = bpp_solver_sgl(ldhw.cols(oinit, otail), m);

      k.subvec(oinit, otail) = afit.k; kd.subvec(oinit, otail) = afit.kb; it.cols(oinit, otail) = afit.it; ok = (ok && afit.ok);

      if (i < id.size()) { oinit = i; olock = id(i); }

    }

  }

  return bppSgl(id, it, bn, k, kd, ok);

}


// bpp_solver_sgl
// @description
//  subroutine of bpp_solver_dpp
// @details
//  fit a single order into bn list
// @return bppSgl
// @family bpp_solver_dpp
// @export
bppSgl bpp_solver_sgl(
  const arma::mat& ldhw, const arma::mat& m
) {

  // init
  arma::uvec id(ldhw.n_cols, arma::fill::zeros); // id - impose id = 0 since packing a single order

  // init
  arma::mat bn = m; // bn

  // init
  arma::mat it(8, ldhw.n_cols, arma::fill::zeros); // it - x, y, z, w (w hold in bn when fit it), l, d, h, w (w of it) on row 0, 1, 2, 3, 4, 5, 6, 7

  arma::uvec xyzw_ulmt = arma::linspace<arma::uvec>(0, 3, 4);

  arma::uvec ldhw_ulmt = arma::linspace<arma::uvec>(4, 7, 4);

  it.rows(xyzw_ulmt).fill(-1); it.rows(ldhw_ulmt) = ldhw; // it - init x, y, z, w, l, d, h, w into it

  // init
  arma::uvec k(ldhw.n_cols, arma::fill::zeros); // k - ticket id vector

  // init
  arma::uvec kb(ldhw.n_cols, arma::fill::zeros); // kb - ticket bn id vector

  // init
  arma::uvec itlmt = arma::linspace<arma::uvec>(0, ldhw.n_cols - 1, ldhw.n_cols); // itlmt: it index no fit yet

  // init
  // bpp_solver_sgl_screen handle corner case when a single it by itself cannot fit into any bn?
  bool ok = bpp_solver_sgl_screen(ldhw, m, itlmt); // impose ticket id 0 if a single it by itself cannot fit into any bn

  // main
  arma::uword tid = 1; // tid: ticket id

  gbp4q afit(arma::zeros<arma::vec>(ldhw.n_cols), it, bn, k, arma::zeros<arma::uvec>(bn.n_cols), 0, ok); arma::uvec klmt; arma::uvec flmt;

  while (itlmt.size() > 0) {

    afit = gbp4d_solver_dpp_filt(ldhw.cols(itlmt), m);

    klmt = arma::find(afit.k == 1); flmt = arma::find(afit.f == 1, 1, "first"); // flmt should always have size 1

    it.cols(itlmt) = afit.it; k(itlmt(klmt)).fill(tid); if (flmt.size() > 0) { kb(itlmt(klmt)).fill(flmt(0)); }

    itlmt = itlmt(arma::find(afit.k != 1)); tid++;

  }

  return bppSgl(id, it, bn, k, kb, ok);
}


// //' bpp_solver_sgl_screen
// //' @description
// //'  auxilium of bpp_solver_sgl
// //' @details
// //'  screen and remove it from itlmt if all(any(sorted (l, d, h, w) > sorted (ml, md, mh, mw))) is true
// //'   or all(v > mv) is true
// //' @return ok?
// //' @family bpp_solver_sgl
// //' @note
// //'  itlmt is in place modified via bpp_solver_sgl_screen
// //' @export
// // [[Rcpp::export]]
bool bpp_solver_sgl_screen(
  const arma::mat& ldhw, const arma::mat& m, arma::uvec& itlmt
) {

  // init
  arma::mat sldhw(4, ldhw.n_cols); // const arma::mat sldhw(4, ldhw.n_cols);

  sldhw.rows(0, 2) = arma::sort(ldhw.rows(0, 2), "descend", 0); sldhw.row(3) = ldhw.row(3);

  arma::mat vldhw(2, ldhw.n_cols); // const arma::vec vldhw(2);

  vldhw.row(0) = ldhw.row(0) % ldhw.row(1) % ldhw.row(2); vldhw.row(1) = ldhw.row(3);

  // m: assume m is sorted l >= d >= h?
  arma::mat sm(4, m.n_cols); // const arma::mat sm(4, m.n_cols);

  sm.rows(0, 2) = arma::sort(m.rows(0, 2), "descend", 0); sm.row(3) = m.row(3);

  arma::mat vm(2, m.n_cols); // const arma::mat vm(2, m.n_cols);

  vm.row(0) = m.row(0) % m.row(1) % m.row(2); vm.row(1) = m.row(3);

  // main
  bool ok = true;

  // an implement assume the last one bn is dominant, check this assumption in bpp_solver_dpp().
  // for (arma::uword i = 0; i < ldh.n_cols; i++) {
  //   if (sldhw(0, i) > sm(0, sm.n_cols - 1) ||
  //       sldhw(1, i) > sm(1, sm.n_cols - 1) ||
  //       sldhw(2, i) > sm(2, sm.n_cols - 1) ||
  //       sldhw(3, i) > sm(3, sm.n_cols - 1) ||
  //       vldhw(0, i) > vm(0, vm.size() - 1) ||
  //       vldhw(1, i) > vm(1, vm.size() - 1) ||
  //   ) {
  //     ok = false; itlmt = itlmt(arma::find(itlmt != i));
  //   }
  // }

  // correct logic however when last bn is not dominant all others can cause forever while loop
  // imagine one scenario when it = c(14, 14, 14), bn0 = c(20, 20, 20), and bn1 = c(50, 20, 10)
  // so we determined it can fit into some bn in bpp_solver_sgl_screen, which actually correct!
  // however, after get back to bpp_solver_sgl() call and into while loop gbp4d_solver_dpp_filt
  // gbp4d_solver_dpp_filt can never fit it into bn because gbp4d_solver_dpp_filt which follows
  // gbp3d_solver_dpp_filt assume last bn is dominant, and it cannot fit into bn, so ok = false
  // and k = 0, and as a result itlmt can never be driven down to zero.

  for (arma::uword i = 0; i < ldhw.n_cols; i++) {
    if (arma::all(sldhw(0, i) > sm.row(0)) ||
        arma::all(sldhw(1, i) > sm.row(1)) ||
        arma::all(sldhw(2, i) > sm.row(2)) ||
        arma::all(sldhw(3, i) > sm.row(3)) ||
        arma::all(vldhw(0, i) > vm.row(0)) // ||
        // arma::all(vldhw(1, i) > vm.row(1)) // same as arma::all(sldhw(3, i) > sm.row(3))
    ) {
      ok = false; itlmt = itlmt(arma::find(itlmt != i));
    }
  }

  return ok;
}
