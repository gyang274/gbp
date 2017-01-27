#include "gbp_u.h"
#include "gbp1d.h"
#include "gbp4d_xp.h"
#include "gbp4d_it.h"
#include "gbp4d.h"


const double tol = 0.00000001;


// class gbp4d {
//
// public:
//
//   arma::vec p; arma::mat it; arma::vec bn; arma::uvec k; double o; bool ok;
//
//   // .constructor
//   gbp4d(
//     arma::vec p, arma::mat it, arma::vec bn, arma::uvec k, double o, bool ok
//   ) : p(p), it(it), bn(bn), k(k), o(o), ok(ok) {
//   }
//
// };

// gbp4d global setup on parameter
const arma::uword gbp4d_nlvl_mkt0 =  2; // when <  2 items search over all in ktlist w.r.t orientation + extreme point
const arma::uword gbp4d_nlvl_mkt1 =  3; // when <  3 items search best   5 in ktlist w.r.t orientation + extreme point
const arma::uword gbp4d_nlvl_mkt2 =  5; // when <  5 items search best   3 in ktlist w.r.t orientation + extreme point
const arma::uword gbp4d_nlvl_mkt3 =  8; // when <  8 items search best   2 in ktlist w.r.t orientation + extreme point
//const arma::uword gbp4d_nlvl_mkt4 = ; // when >= 8 items search best   1 in ktlist w.r.t orientation + extreme point

// solve gbp4d via extreme point heuristic and best information score fit strategy
// gbp4d gbp4d_solver_dpp(const arma::vec& p, const arma::mat& ldhw, const arma::vec& m);
gbp4d gbp4d_solver_dpp(const arma::vec& p, const arma::mat& ldhw, const arma::vec& m) {

  // init
  // void gbp4d_solver_dpp_init() {}

  // init parameters in global solver
  const arma::uword n = p.size(); // n: number of it

  // if (n != ldhw.n_cols) { throw std::exception(); }

  const arma::uvec q = arma::sort_index(p, "descend"); // q: it fit sequence

  const arma::vec bn = m; // bn

  arma::mat it = arma::zeros<arma::mat>(8, n); // it
  arma::uvec xyzw_ulmt = arma::linspace<arma::uvec>(0, 3, 4); // it - x, y, z, w on row 0, 1, 2, 3
  arma::uvec ldhw_ulmt = arma::linspace<arma::uvec>(4, 7, 4); // it - l, d, h, w on row 4, 5, 6, 7
  it.rows(ldhw_ulmt) = ldhw; // it - init l, d, h into it

  arma::mat itastr = it; // itastr: it asterisk: global track current best it over recursive call

  arma::mat xp = arma::zeros<arma::mat>(8, 1); xp(4, 0) = bn(0); xp(5, 0) = bn(1); xp(6, 0) = bn(2); xp(7, 0) = bn(3); // xp

  arma::uword nlvl = n; // nlvl: number of it open to fit

  arma::uword nastr = gbp4d_solver_dpp_main_create_nastr(p, ldhw, m); // nastr: number of it oversize volume or weight limit

  arma::uvec g = arma::zeros<arma::uvec>(n); // g: status vector: 0 open, 1 fitted, and 2 determined no fit

  arma::uvec gastr = g; // gastr: g asterisk: global track current best g over recursive call

  const arma::vec v = arma::trans(it.row(4) % it.row(5) % it.row(6)); // v: volume of it - l d h in objective and constraint while w in constraint only

  const double vastr = std::min(arma::sum(v), bn(0) * bn(1) * bn(2)); // vastr: v asterisk: maximum volume can be achieved via fit

  double u = 0; // u: sum of v volume fitted into bn

  double uastr = u; // uastr: u asterisk: global track current best u over recursive call

  arma::uvec k = arma::zeros<arma::uvec>(n); // k: selection indicator <vector>

  double o = 0; // o: objective achievement <double>

  // main fit corner case when n == 0
  if (n == 0 || m.n_cols == 0) { return gbp4d(p, it, bn, k, o, true); }

  // main fit recursive call - should create a class gbp4d_allinfo hold all info?
  bool ok = gbp4d_solver_dpp_main(bn, it, itastr, xp, q, nlvl, nastr, g, gastr, v, vastr, u, uastr);

  // ok via gbp4d_solver_dpp(): can all it fit into bn?
  // ok via gbp4d_solver_dpp_main() hold a different meaning for recursive purpose
  // ok via gbp4d_solver_dpp_main(): can all it fit into bn or a subset of it make full use of bn 100% volume utilization?
  arma::uvec glmt = arma::find(gastr == 1); k(glmt).fill(1); o = arma::sum(v(glmt)); ok = arma::all(gastr == 1);

  // it via itastr and gastr: flag no fit with (x, y) = (-1, -1, -1) instead of (0, 0, 0)
  it.rows(xyzw_ulmt).fill(-1); it.cols(glmt) = itastr.cols(glmt);

  return gbp4d(p, it, bn, k, o, ok);
}

// //' gbp4d_solver_dpp_main
// //' @description
// //'  gbp4d_solver_dpp main recursive solver
// //' @details
// //'  create it fit via recursive fit it into bn
// //' @return ok
// //' @family gbp4d_solver_dpp
// //' @export
// // [[Rcpp::export]]
bool gbp4d_solver_dpp_main(
    const arma::vec& bn, arma::mat& it, arma::mat& itastr, arma::mat& xp, const arma::uvec& q, const arma::uword nlvl, const arma::uword nastr,
    arma::uvec g, arma::uvec& gastr, const arma::vec& v, const double& vastr, double u, double& uastr
) {

  // init

  arma::uword id = q(q.size() - nlvl); // index of it to fit

  // index of g into ok - all it not fit or determined not fit yet
  arma::uvec glmt = arma::zeros<arma::uvec>(g.size()); glmt(arma::find(g == 0)) += 1;

  bool ok = false; // ok - can all it fit into bn or a subset of it make full use of bn 100% volume utilization

  // create ktlist
  arma::mat it0 = it.cols(arma::find(g == 1)); arma::mat xp0 = xp; arma::vec kt0 = it.col(id);

  arma::uword nlmt = gbp4d_solver_dpp_main_create_nlmt(nlvl, nastr);

  Ktlist4d ktlist = gbp4d_it_create_ktlist(bn, it0, xp0, kt0, nlmt);

  // main
  if (nlvl == 1) {

    // fit final it in queue
    if (ktlist.n > 0) {

      it.col(id) = ktlist.kt.col(0); g(id) = 1; u = u + v(id);

      if (u > uastr) {

        itastr = it; gastr = g; uastr = u;

      }

      if (uastr == vastr || arma::all(gastr == 1)) {

        ok = true; return ok;

      }

    } else {

      g(id) = 2;

    }

  } else {

    // do recursive call
    if (ktlist.n > 0) {

      g(id) = 1; u = u + v(id);

      for (arma::uword i = 0; i < ktlist.n; i++) {

        it.col(id) = ktlist.kt.col(i); xp = ktlist.xp(i);

        if (u > uastr) {

          itastr = it; gastr = g; uastr = u;

        }

        if (uastr == vastr || arma::all(gastr == 1)) {

          ok = true; return ok;

        }

        ok = gbp4d_solver_dpp_main(bn, it, itastr, xp, q, nlvl - 1, nastr, g, gastr, v, vastr, u, uastr);

        if (ok) { break; }

      }

      // what if skip this it? can later it combined be better?
      if (!ok && nlvl < gbp4d_nlvl_mkt2) {

        double vmiss = arma::sum(v(arma::find(gastr != 1 and glmt == 1)));

        if (vmiss > v(id)) {

          g(id) = 2; u = u - v(id); it(0, id) = 0; it(1, id) = 0; it(2, id) = 0; it(3, id) = 0; xp = xp0;

          ok = gbp4d_solver_dpp_main(bn, it, itastr, xp, q, nlvl - 1, nastr, g, gastr, v, vastr, u, uastr);

        }

      }

    } else {

      g(id) = 2; ok = gbp4d_solver_dpp_main(bn, it, itastr, xp, q, nlvl - 1, nastr, g, gastr, v, vastr, u, uastr);

    }

  }

  return ok;
}

// //' gbp4d_solver_dpp_main_create_nastr
// //' @description
// //'  subroutine of gbp4d_solver_dpp_main
// //' @details
// //'  create nastr
// //' @return nastr
// //' @family gbp4d_solver_dpp
// //' @export
// // [[Rcpp::export]]
arma::uword gbp4d_solver_dpp_main_create_nastr(const arma::vec& p, const arma::mat& ldhw, const arma::vec& m) {

  // init
  const arma::uvec q = arma::sort_index(p, "descend");

  const arma::vec v = arma::trans(ldhw.row(0) % ldhw.row(1) % ldhw.row(2));

  const arma::vec w = arma::trans(ldhw.row(3));

  const double mv = arma::prod(m.subvec(0, 2));

  const double mw = m(3);

  // main
  double v0 = 0.00; double w0 = 0.00;

  arma::uword nastr = 0;

  for (arma::uword i = 0; i < q.size(); i++) {

    v0 += v(q(i)); w0 += w(q(i));

    if ((v0 >= mv) || (w0 >= mw)) {

      nastr = q.size() - 1 - i; break;

    }

  }

  return nastr;
}

// //' gbp4d_solver_dpp_main_create_nlmt
// //' @description
// //'  subroutine of gbp4d_solver_dpp_main
// //' @details
// //'  create nlmt in ktlist
// //' @return nlmt
// //' @family gbp4d_solver_dpp
// //' @export
// // [[Rcpp::export]]
arma::uword gbp4d_solver_dpp_main_create_nlmt(const arma::uword nlvl, const arma::uword nastr) {

  arma::uword nlmt = 1;

  if (nlvl <= nastr) { return nlmt; }

  if (nlvl - nastr < gbp4d_nlvl_mkt0) {
    nlmt = 0;
  } else if (nlvl - nastr < gbp4d_nlvl_mkt1) {
    nlmt = 5;
  } else if (nlvl - nastr < gbp4d_nlvl_mkt2) {
    nlmt = 3;
  } else if (nlvl - nastr < gbp4d_nlvl_mkt3) {
    nlmt = 2;
  } else {
    nlmt = 1;
  }

  return nlmt;
}

//' gbp4d_solver_dpp_prep_create_p
//' @description
//'  auxilium of gbp4d_solver_dpp
//' @details
//'  create p via ldhw and m via cluster w, cluster max(l, d, h) and area strategy
//' @param ldhw 4xN matrix of l, d, h, w of it
//' @param m 4x1 vector of l, d, h, w of bn
//' @return p
//' @family gbp4d_solver_dpp
//' @export
// [[Rcpp::export]]
arma::vec gbp4d_solver_dpp_prep_create_p(const arma::mat& ldhw, const arma::vec& m) {

  // init
  arma::vec p = arma::zeros<arma::vec>(ldhw.n_cols);

  if (ldhw.n_cols == 0 || m.size() == 0) { return p; }

  // main
  arma::mat q(3, ldhw.n_cols);

  // create it w cluster w.r.t gbp1d
  arma::vec v = arma::trans(ldhw.row(0) % ldhw.row(1) % ldhw.row(2));

  arma::uvec w = arma::conv_to<arma::uvec>::from(arma::trans(arma::ceil(ldhw.row(3) / 0.25)));

  arma::uword c = std::floor(m(3) / 0.25);

  gbp1d wfit = gbp1d_solver_dpp(v, w, c);

  for (arma::uword i = 0; i < ldhw.n_cols; i++) {

    if (wfit.k(i) == 1) {

      q(0, i) = 1.0;

    } else {

      q(0, i) = 0.0;

    }

  }

  // create it max(l, d, h) cluster w.r.t bn max(l, d, h) into per 0.25 inch
  double mldh = m.subvec(0, 2).max(); arma::rowvec ildh = arma::max(ldhw.rows(0, 2), 0);

  arma::vec nldh = arma::linspace<arma::vec>(0, mldh, std::ceil(mldh * 4.0));

  for (arma::uword i = 0; i < ldhw.n_cols; i++) {

    q(1, i) = arma::sum(arma::find(nldh <= ildh(i)));

  }

  q.row(2) = ldhw.row(0) % ldhw.row(1) % ldhw.row(2) / ildh;

  // sort index after sort index via rows so that it with higher value of p can get fit into bn earlier in queue
  p = arma::conv_to<arma::vec>::from(sort_index(sort_index_via_rows(q, arma::linspace<arma::uvec>(0, 2, 3))));

  return p;
}


// class gbp4q {
//
// public:
//
//   arma::vec p; arma::mat it; arma::mat bn; arma::uvec k; arma::uvec f; double o; bool ok;
//
//   // .constructor
//   gbp4q(
//     arma::vec p, arma::mat it, arma::mat bn, arma::uvec k, arma::uvec f, double o, bool ok
//   ) : p(p), it(it), bn(bn), k(k), f(f), o(o), ok(ok) {
//   }
//
// };

// gbp4d_solver_dpp_filt
// @description
//  auxilium of gbp4d_solver_dpp
// @param ldhw
//  it scale <matrix>
//  - l, d, h, w it scale along x, y, z and also w <numeric>
// @param m
//  bn scale <matrix>
//  - l, d, h, w bn scale along x, y, z and also w <numeric>
//  - l, d, h, w in row and each col is a single bn
//  should make sure bn list are sorted via volume
//   so that the first col is the most prefered smallest bn, and also
//   the last col is the least prefered largest and often dominant bn
//  should make sure no X in front of Y if bnX dominant bnY,
//   bnX dominant bnY if all(X(l, d, h) > Y(l, d, h)) and should always prefer Y.
//  should make sure bn such that l >= d >= h or vice versa.
// @details
//  gbp4d_solver_dpp_filt is built on top of gbp4d_solver_dpp
//   aims to select the most preferable bn from a list of bn that can fit all or most it
//  gbp4d_solver_dpp()'s objective is fit all or most it into a single given bn (l, d, h)
//  gbp4d_solver_dpp_filt()'s objective is select the most preferable given a list of bn
//   where bn list is specified in 3xN matrix that the earlier column the more preferable
//  gbp4d_solver_dpp_filt() use an approx binary search and determine f w.r.t bn.n_cols
//   where f = 1 indicate the bn being selected and only one of 1 in result returned.
//  ok = true if any bin can fit all it and algorithm will select smallest bn can fit all
//   otherwise ok = false and algorithm will select a bn can maximize volume of fitted it
//  often recommend to make the last and least preferable bn dominate all other bn in list
//   when design bn list, bnX dominant bnY if all(X(l, d, h) > Y(l, d, h)).
// @return gbp4q
// @family gbp4d_solver_dpp_filt
// @export
gbp4q gbp4d_solver_dpp_filt(const arma::mat& ldhw, const arma::mat& m) {

  // init
  arma::mat sldhw(4, ldhw.n_cols); // const arma::mat sldhw(4, ldhw.n_cols);

  sldhw.rows(0, 2) = arma::sort(ldhw.rows(0, 2), "descend", 0); sldhw.row(3) = ldhw.row(3);

  arma::vec vldhw(2); // const arma::vec vldhw(2);

  vldhw(0) = arma::sum(ldhw.row(0) % ldhw.row(1) % ldhw.row(2)); vldhw(1) = arma::sum(ldhw.row(3));

  // m: assume m is sorted l >= d >= h?
  arma::mat sm(4, m.n_cols); // const arma::mat sm(4, m.n_cols);

  sm.rows(0, 2) = arma::sort(m.rows(0, 2), "descend", 0); sm.row(3) = m.row(3);

  arma::mat vm(2, m.n_cols); // const arma::mat vm(2, m.n_cols);

  vm.row(0) = m.row(0) % m.row(1) % m.row(2); vm.row(1) = m.row(3);

  // main
  arma::uvec flmt = arma::linspace<arma::uvec>(0, m.n_cols - 1, m.n_cols);

  // main: assume ldhw.n_cols > 0 and m.n_cols > 0

  // main: corner case: ldhw.n_cols == 1 fast filt
  if (ldhw.n_cols == 1) {

    gbp4d_solver_dpp_filt_fast(sldhw, vldhw, sm, vm, flmt);

  }

  // main: id p fit: global parameter track current selected bn from bn list;
  // main: init fit it into last bin if last largest bin that can often drive down search space
  // ldh itlmt should be pre-screened via bpp_solver_sgl_screen so that flmt.size() always > 0 at this stage
  arma::uword id = m.n_cols - 1; if (flmt.size() > 0) { id = flmt(flmt.size() - 1); flmt = flmt(arma::find(flmt != id)); }

  arma::vec p = gbp4d_solver_dpp_prep_create_p(ldhw, m.col(id)); gbp4d fit = gbp4d_solver_dpp(p, ldhw, m.col(id));

  if (fit.ok) {

    gbp4d_solver_dpp_filt_fast(sldhw, vldhw, sm, vm, flmt); // fast filt g given known all it can fit into some one bn

  } else {

    gbp4d_solver_dpp_filt_slow(fit.o, vm.row(0), flmt); // slow filt g given known not all but at least some it can fit into some one bn

  }

  // main: fit it into each single bn candidate via approx binary search
  arma::uword id0; arma::vec p0; gbp4d fit0(arma::vec(), arma::mat(), arma::vec(), arma::uvec(), 0.00, false);

  while(flmt.size() > 0) {

    // id0: must use floor since in corner case flmt.size() == 1 can only take id0 = flmt(0);
    id0 = flmt(std::floor(flmt.size() / 2.00));

    p0 = gbp4d_solver_dpp_prep_create_p(ldhw, m.col(id0)); fit0 = gbp4d_solver_dpp(p0, ldhw, m.col(id0));

    if (fit0.ok) {

      // run gbp4d_solver_dpp_filt_fast once when fit.ok becomes ok == true;
      if (!fit.ok) {

        gbp4d_solver_dpp_filt_fast(sldhw, vldhw, sm, vm, flmt);

      }

      // assume m column is sorted via volume, otherwise vm(flmt) > vm(id0)
      flmt = flmt(arma::find(flmt < id0)); id = id0; p = p0; fit = fit0; // if (vm(id0) <= vm(id)) { id = id0; p = p0; fit = fit0; }

    } else {

      gbp4d_solver_dpp_filt_slow(fit.o, vm.row(0), flmt);

      if (fit.ok) {

        gbp4d_solver_dpp_filt_dcol(sm.rows(0, 2), id0, flmt);

      } else {

        // assume m column is sorted via volume, otherwise vm(id0) < vm(id)
        if (fit0.o > fit.o || (fit0.o == fit.o && id0 < id)) { id = id0; p = p0; fit = fit0; }

        // corner case protection
        // it0 - it9, bn0 and bn1, and bn0 is smaller than bn1
        // fit it0 - it8 into bn1, fit0 it0 - it7 and it9 into bn0, sum(v(it0 - it7, it9)) < sum(v(it0 - it8)),
        // while it is possible to also fit it0 - it8 into bn0 with recursive, it didn't get a chance recursive
        // since it9 is available and just fit, so a reinforcement recursive search if fit0 close enough to fit
        // if (fit0.o < fit.o && id0 < id) {
        //
        //   arma::uvec klmt = arma::find(fit.k == 1); fit0 = gbp4d_solver_dpp(p0(klmt), ldhw.cols(klmt), m.col(id0));
        //
        //   if (fit0.ok) { id = id0; p = p0; fit.it.cols(klmt) = fit0.it; }
        //
        // }

      }

    }

    flmt = flmt(arma::find(flmt != id0));

  }

  arma::uvec f(m.n_cols, arma::fill::zeros); f(id) = 1;

  return gbp4q(p, fit.it, m, fit.k, f, fit.o, fit.ok);
}

// //' assert_matrix_last_col_dominant
// //' @description
// //'  auxilium of gbp4d_solver_dpp_filt
// //' @details
// //'  assert matrix last column dominant all other column: arma::all(m.col(i) <= m.col(m.n_cols - 1)), all i.
// //' @return ok_last_col_dominant
// //' @family gbp4d_solver_dpp_filt
// //' @export
// // [[Rcpp::export]]
// bool assert_matrix_last_col_dominant(const arma::mat& m) {
//
//   bool ok_last_col_dominant = true;
//
//   for (arma::uword i = 0; i < m.n_cols - 1; i++) {
//     // should make sure m such that l >= d >= h or vice versa and w on separate single dimension
//     if (arma::all(m.col(i) <= m.col(m.n_cols - 1) + tol)) { ok_last_col_dominant = false; break; }
//   }
//
//   return ok_last_col_dominant;
// }

// //' gbp4d_solver_dpp_filt_fast
// //' @description
// //'  auxilium of gbp4d_solver_dpp_filt
// //' @details
// //'  gbp4d_solver_dpp_filt_fast use in context of gbp4d_solver_dpp_filt:
// //'   when determine all it can fit into some one bn, then can safely filt
// //'   out bn such that cannot hold all it at once determined via comparing
// //'   scale, volume and weight
// //' @return flmt
// //' @family gbp4d_solver_dpp_filt
// //' @note
// //'  call in r won't see result since flmt is changed in place - designed to use in rcpp
// //' @export
// // [[Rcpp::export]]
void gbp4d_solver_dpp_filt_fast(
    const arma::mat& sldhw, const arma::vec& vldhw, const arma::mat& sm, const arma::mat& vm, arma::uvec& flmt
) {

  // main
  arma::uvec flmtidk(flmt.size(), arma::fill::zeros);

  for (arma::uword i = 0; i < flmt.size(); i++) {
    if (arma::any(sldhw.row(0) > sm(0, flmt(i))) ||
        arma::any(sldhw.row(1) > sm(1, flmt(i))) ||
        arma::any(sldhw.row(2) > sm(2, flmt(i))) ||
        arma::any(sldhw.row(3) > sm(3, flmt(i))) ||
        vldhw(0) > vm(0, flmt(i)) ||
        vldhw(1) > vm(1, flmt(i))
    ) {
      flmtidk(i) = 1;
    }
  }

  flmt = flmt(arma::find(flmtidk != 1));

}

// //' gbp4d_solver_dpp_filt_slow
// //' @description
// //'  auxilium of gbp4d_solver_dpp_filt
// //' @details
// //' gbp4d_solver_dpp_filt_slow use in context of gbp4d_solver_dpp_filt:
// //'  when determine some it with volume v0 can fit into some one bn,
// //'  can safely filt out bn such that the bn volume is less than v0.
// //' @return flmt
// //' @family gbp4d_solver_dpp_filt
// //' @note
// //'  call in r won't see result since flmt is changed in place - designed to use in rcpp
// //' @export
// // [[Rcpp::export]]
void gbp4d_solver_dpp_filt_slow(
    const double vldh, const arma::rowvec& vm, arma::uvec& flmt
) {

  // main
  flmt = flmt(arma::find(vm(flmt) >= vldh)); // assume bn w is none decreasing with bn v = l * d * h since w is not considered as objective but constraint only

}

// //' gbp4d_solver_dpp_filt_dcol
// //' @description
// //'  auxilium of gbp4d_solver_dpp_filt
// //' @details
// //'  gbp4d_solver_dpp_filt_dcol use in context of gbp4d_solver_dpp_filt:
// //'   when determine all it can fit into some one bn, and a particular bn bn0
// //'   cannot fit all it, then can safely filt out bn1 dominated by bn0 in the
// //'   sense that bn1 also cannot fit all it.
// //'  bn1 dominated by bn0 if all bn1(l, d, h) <= bn0(l, d, h), l >= d >= h.
// //'  must first determine all it can fit into some one bn, otherwise it is possible
// //'   that no bn can fit all it, so even it is correct that if bn0 cannot fit all,
// //'   then bn1 cannot neither, but bn1 is smaller than bn0, and could fit the same
// //'   set of it and thus provide a higher preference w.r.t utilization rate.
// //' @return flmt
// //' @family gbp4d_solver_dpp_filt
// //' @note
// //'  call in r won't see result since g is changed in place - designed to use in rcpp
// //' @export
// // [[Rcpp::export]]
void gbp4d_solver_dpp_filt_dcol(
    const arma::mat &sm, const arma::uword icol, arma::uvec& flmt
) {

  // main
  arma::uvec flmtidk(flmt.size(), arma::fill::zeros);

  for (arma::uword i = 0; i < flmt.size(); i++) {
    if (arma::all(sm.col(flmt(i)) <= sm.col(icol) + tol)) {
      flmtidk(i) = 1;
    }
  }

  flmt = flmt(arma::find(flmtidk != 1));

}
