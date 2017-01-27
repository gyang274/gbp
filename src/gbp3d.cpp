#include "gbp_u.h"
#include "gbp3d_xp.h"
#include "gbp3d_it.h"
#include "gbp3d.h"


const double tol = 0.00000001;


// class gbp3d {
//
// public:
//
//   arma::vec p; arma::mat it; arma::vec bn; arma::uvec k; double o; bool ok;
//
//   // .constructor
//   gbp3d(
//     arma::vec p, arma::mat it, arma::vec bn, arma::uvec k, double o, bool ok
//   ) : p(p), it(it), bn(bn), k(k), o(o), ok(ok) {
//   }
//
// };

// gbp3d global setup on parameter
const arma::uword gbp3d_nlvl_mkt0 =  2; // when <  2 items search over all in ktlist w.r.t orientation + extreme point
const arma::uword gbp3d_nlvl_mkt1 =  3; // when <  3 items search best   5 in ktlist w.r.t orientation + extreme point
const arma::uword gbp3d_nlvl_mkt2 =  5; // when <  5 items search best   3 in ktlist w.r.t orientation + extreme point
const arma::uword gbp3d_nlvl_mkt3 =  8; // when <  8 items search best   2 in ktlist w.r.t orientation + extreme point
//const arma::uword gbp3d_nlvl_mkt4 = ; // when >= 8 items search best   1 in ktlist w.r.t orientation + extreme point

// solve gbp3d via extreme point heuristic and best information score fit strategy
// gbp3d gbp3d_solver_dpp(const arma::vec& p, const arma::mat& ldh, const arma::vec& m);
gbp3d gbp3d_solver_dpp(const arma::vec& p, const arma::mat& ldh, const arma::vec& m) {

  // init
  // void gbp3d_solver_dpp_init() {}

  // init parameters in global solver
  const arma::uword n = p.size(); // n: number of it

  // if (n != ldh.n_cols) { throw std::exception(); }

  const arma::uvec q = arma::sort_index(p, "descend"); // q: it fit sequence

  const arma::vec bn = m; // bn

  arma::mat it = arma::zeros<arma::mat>(6, n); // it
  arma::uvec xyz_ulmt = arma::linspace<arma::uvec>(0, 2, 3); // it - x, y, z on row 0, 1, 2
  arma::uvec ldh_ulmt = arma::linspace<arma::uvec>(3, 5, 3); // it - l, d, h on row 3, 4, 5
  it.rows(ldh_ulmt) = ldh; // it - init l, d, h into it

  arma::mat itastr = it; // itastr: it asterisk: global track current best it over recursive call

  arma::mat xp = arma::zeros<arma::mat>(6, 1); xp(3, 0) = bn(0); xp(4, 0) = bn(1); xp(5, 0) = bn(2); // xp

  arma::uword nlvl = n; // nlvl: number of it open to fit

  arma::uword nastr = gbp3d_solver_dpp_main_create_nastr(p, ldh, m); // nastr: number of it oversize volume or weight limit

  arma::uvec g = arma::zeros<arma::uvec>(n); // g: status vector: 0 open, 1 fitted, and 2 determined no fit

  arma::uvec gastr = g; // gastr: g asterisk: global track current best g over recursive call

  const arma::vec v = arma::trans(it.row(3) % it.row(4) % it.row(5)); // v: volume of it

  const double vastr = std::min(arma::sum(v), bn(0) * bn(1) * bn(2)); // vastr: v asterisk: maximum volume can be achieved via fit

  double u = 0; // u: sum of v volume fitted into bn

  double uastr = u; // uastr: u asterisk: global track current best u over recursive call

  arma::uvec k = arma::zeros<arma::uvec>(n); // k: selection indicator <vector>

  double o = 0; // o: objective achievement <double>

  // main fit corner case when n == 0
  if (n == 0 || m.n_cols == 0) { return gbp3d(p, it, bn, k, o, true); }

  // main fit recursive call - should create a class gbp3d_allinfo hold all info?
  bool ok = gbp3d_solver_dpp_main(bn, it, itastr, xp, q, nlvl, nastr, g, gastr, v, vastr, u, uastr);

  // ok via gbp3d_solver_dpp(): can all it fit into bn?
  // ok via gbp3d_solver_dpp_main() hold a different meaning for recursive purpose
  // ok via gbp3d_solver_dpp_main(): can all it fit into bn or a subset of it make full use of bn 100% volume utilization?
  arma::uvec glmt = arma::find(gastr == 1); k(glmt).fill(1); o = arma::sum(v(glmt)); ok = arma::all(gastr == 1);

  // it via itastr and gastr: flag no fit with (x, y) = (-1, -1, -1) instead of (0, 0, 0)
  it.rows(xyz_ulmt).fill(-1); it.cols(glmt) = itastr.cols(glmt);

  return gbp3d(p, it, bn, k, o, ok);
}

// //' gbp3d_solver_dpp_main
// //' @description
// //'  gbp3d_solver_dpp main recursive solver
// //' @details
// //'  create it fit via recursive fit it into bn
// //' @return ok
// //' @family gbp3d_solver_dpp
// //' @export
// // [[Rcpp::export]]
bool gbp3d_solver_dpp_main(
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

  arma::uword nlmt = gbp3d_solver_dpp_main_create_nlmt(nlvl, nastr);

  Ktlist3d ktlist = gbp3d_it_create_ktlist(bn, it0, xp0, kt0, nlmt);

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

        ok = gbp3d_solver_dpp_main(bn, it, itastr, xp, q, nlvl - 1, nastr, g, gastr, v, vastr, u, uastr);

        if (ok) { break; }

      }

      // what if skip this it? can later it combined be better?
      if (!ok && nlvl < gbp3d_nlvl_mkt2) {

        double vmiss = arma::sum(v(arma::find(gastr != 1 and glmt == 1)));

        if (vmiss > v(id)) {

          g(id) = 2; u = u - v(id); it(0, id) = 0; it(1, id) = 0; it(2, id) = 0; xp = xp0;

          ok = gbp3d_solver_dpp_main(bn, it, itastr, xp, q, nlvl - 1, nastr, g, gastr, v, vastr, u, uastr);

        }

      }

    } else {

      g(id) = 2; ok = gbp3d_solver_dpp_main(bn, it, itastr, xp, q, nlvl - 1, nastr, g, gastr, v, vastr, u, uastr);

    }

  }

  return ok;
}

// //' gbp3d_solver_dpp_main_create_nastr
// //' @description
// //'  subroutine of gbp3d_solver_dpp_main
// //' @details
// //'  create nastr
// //' @return nastr
// //' @family gbp3d_solver_dpp
// //' @export
// // [[Rcpp::export]]
arma::uword gbp3d_solver_dpp_main_create_nastr(const arma::vec& p, const arma::mat& ldh, const arma::vec& m) {

  // init
  const arma::uvec q = arma::sort_index(p, "descend");

  const arma::vec v = arma::trans(ldh.row(0) % ldh.row(1) % ldh.row(2));

  const double mv = arma::prod(m);

  // main
  double v0 = 0.00;

  arma::uword nastr = 0;

  for (arma::uword i = 0; i < q.size(); i++) {

    v0 += v(q(i));

    if ((v0 >= mv)) {

      nastr = q.size() - 1 - i; break;

    }

  }

  return nastr;
}

// //' gbp3d_solver_dpp_main_create_nlmt
// //' @description
// //'  subroutine of gbp3d_solver_dpp_main
// //' @details
// //'  create nlmt in ktlist
// //' @return nlmt
// //' @family gbp3d_solver_dpp
// //' @export
// // [[Rcpp::export]]
arma::uword gbp3d_solver_dpp_main_create_nlmt(const arma::uword nlvl, const arma::uword nastr) {

  arma::uword nlmt = 1;

  if (nlvl <= nastr) { return nlmt; }

  if (nlvl - nastr < gbp3d_nlvl_mkt0) {
    nlmt = 0;
  } else if (nlvl - nastr < gbp3d_nlvl_mkt1) {
    nlmt = 5;
  } else if (nlvl - nastr < gbp3d_nlvl_mkt2) {
    nlmt = 3;
  } else if (nlvl - nastr < gbp3d_nlvl_mkt3) {
    nlmt = 2;
  } else {
    nlmt = 1;
  }

  return nlmt;
}

//' gbp3d_solver_dpp_prep_create_p
//' @description
//'  auxilium of gbp3d_solver_dpp
//' @details
//'  create p via ldh and m via cluster max(l, d, h) and area strategy
//' @param ldh 3xN matrix of l, d, h of it
//' @param m 3x1 vector of l, d, h of bn
//' @return p
//' @family gbp3d_solver_dpp
//' @export
// [[Rcpp::export]]
arma::vec gbp3d_solver_dpp_prep_create_p(const arma::mat& ldh, const arma::vec& m) {

  // init
  arma::vec p = arma::zeros<arma::vec>(ldh.n_cols);

  if (ldh.n_cols == 0 || m.size() == 0) { return p; }

  // main
  arma::mat q(2, ldh.n_cols);

  // create it max(l, d, h) cluster w.r.t bn max(l, d, h) into per 0.25 inch
  double mldh = m.max(); arma::rowvec ildh = arma::max(ldh, 0);

  arma::vec nldh = arma::linspace<arma::vec>(0, mldh, std::ceil(mldh * 4.0));

  for (arma::uword i = 0; i < ldh.n_cols; i++) {

    q(0, i) = arma::sum(arma::find(nldh <= ildh(i)));

  }

  q.row(1) = ldh.row(0) % ldh.row(1) % ldh.row(2) / ildh;

  // sort index after sort index via rows so that it with higher value of p can get fit into bn earlier in queue
  p = arma::conv_to<arma::vec>::from(sort_index(sort_index_via_rows(q, arma::linspace<arma::uvec>(0, 1, 2))));

  return p;
}


// class gbp3q {
//
// public:
//
//   arma::vec p; arma::mat it; arma::mat bn; arma::uvec k; arma::uvec f; double o; bool ok;
//
//   // .constructor
//   gbp3q(
//     arma::vec p, arma::mat it, arma::mat bn, arma::uvec k, arma::uvec f, double o, bool ok
//   ) : p(p), it(it), bn(bn), k(k), f(f), o(o), ok(ok) {
//   }
//
// };

// gbp3d_solver_dpp_filt
// @description
//  auxilium of gbp3d_solver_dpp
// @param ldh
//  it scale <matrix>
//  - l, d, h it scale along x, y, z <numeric>
// @param m
//  bn scale <matrix>
//  - l, d, h bn scale along x, y, z <numeric>
//  - l, d, h in row and each col is a single bn
//  should make sure bn list are sorted via volume
//   so that the first col is the most prefered smallest bn, and also
//   the last col is the least prefered largest and often dominant bn
//  should make sure no X in front of Y if bnX dominant bnY,
//   bnX dominant bnY if all(X(l, d, h) > Y(l, d, h)) and should always prefer Y.
//  should make sure bn such that l >= d >= h or vice versa.
// @details
//  gbp3d_solver_dpp_filt is built on top of gbp3d_solver_dpp
//   aims to select the most preferable bn from a list of bn that can fit all or most it
//  gbp3d_solver_dpp()'s objective is fit all or most it into a single given bn (l, d, h)
//  gbp3d_solver_dpp_filt()'s objective is select the most preferable given a list of bn
//   where bn list is specified in 3xN matrix that the earlier column the more preferable
//  gbp3d_solver_dpp_filt() use an approx binary search and determine f w.r.t bn.n_cols
//   where f = 1 indicate the bn being selected and only one of 1 in result returned.
//  ok = true if any bin can fit all it and algorithm will select smallest bn can fit all
//   otherwise ok = false and algorithm will select a bn can maximize volume of fitted it
//  often recommend to make the last and least preferable bn dominate all other bn in list
//   when design bn list, bnX dominant bnY if all(X(l, d, h) > Y(l, d, h)).
// @return gbp3q
// @family gbp3d_solver_dpp_filt
// @export
gbp3q gbp3d_solver_dpp_filt(const arma::mat& ldh, const arma::mat& m) {

  // init
  const arma::mat sldh = arma::sort(ldh, "descend", 0);

  const double vldh = arma::sum(ldh.row(0) % ldh.row(1) % ldh.row(2));

  // m: assume m is sorted l >= d >= h?
  const arma::mat sm = arma::sort(m, "descend", 0);

  const arma::rowvec vm = m.row(0) % m.row(1) % m.row(2);

  // main
  arma::uvec flmt = arma::linspace<arma::uvec>(0, m.n_cols - 1, m.n_cols);

  // main: assume ldh.n_cols > 0 and m.n_cols > 0

  // main: corner case: ldh.n_cols == 1 fast filt
  if (ldh.n_cols == 1) {

    gbp3d_solver_dpp_filt_fast(sldh, vldh, sm, vm, flmt);

  }

  // main: id p fit: global parameter track current selected bn from bn list;
  // main: init fit it into last bin if last largest bin that can often drive down search space
  // ldh itlmt should be pre-screened via bpp_solver_sgl_screen so that flmt.size() always > 0 at this stage
  arma::uword id = m.n_cols - 1; if (flmt.size() > 0) { id = flmt(flmt.size() - 1); }

  arma::vec p = gbp3d_solver_dpp_prep_create_p(ldh, m.col(id)); gbp3d fit = gbp3d_solver_dpp(p, ldh, m.col(id));

  if (fit.ok) {

    gbp3d_solver_dpp_filt_fast(sldh, vldh, sm, vm, flmt); // fast filt g given known all it can fit into some one bn

  } else {

    gbp3d_solver_dpp_filt_slow(fit.o, vm, flmt); // slow filt g given known not all but at least some it can fit into some one bn

  }

  // main: fit it into each single bn candidate via approx binary search
  arma::uword id0; arma::vec p0; gbp3d fit0(arma::vec(), arma::mat(), arma::vec(), arma::uvec(), 0.00, false);

  while(flmt.size() > 0) {

    // id0: must use floor since in corner case flmt.size() == 1 can only take id0 = flmt(0);
    id0 = flmt(std::floor(flmt.size() / 2.00));

    p0 = gbp3d_solver_dpp_prep_create_p(ldh, m.col(id0)); fit0 = gbp3d_solver_dpp(p0, ldh, m.col(id0));

    if (fit0.ok) {

      // run gbp3d_solver_dpp_filt_fast once when fit.ok becomes ok == true;
      if (!fit.ok) {

        gbp3d_solver_dpp_filt_fast(sldh, vldh, sm, vm, flmt);

      }

      // assume m column is sorted via volume, otherwise vm(flmt) > vm(id0)
      flmt = flmt(arma::find(flmt < id0)); id = id0; p = p0; fit = fit0; // if (vm(id0) <= vm(id)) { id = id0; p = p0; fit = fit0; }

    } else {

      gbp3d_solver_dpp_filt_slow(fit.o, vm, flmt);

      if (fit.ok) {

        gbp3d_solver_dpp_filt_dcol(sm, id0, flmt);

      } else {

        // assume m column is sorted via volume, otherwise vm(id0) < vm(id)
        if (fit0.o > fit.o || (fit0.o == fit.o && id0 < id)) { id = id0; p = p0; fit = fit0; }

        // if (fit0.o < fit.o && id0 < id) {
        //
        //   arma::uvec klmt = arma::find(fit.k == 1); fit0 = gbp3d_solver_dpp(p0(klmt), ldh.cols(klmt), m.col(id0));
        //
        //   if (fit0.ok) { id = id0; p = p0; fit.it.cols(klmt) = fit0.it; }
        //
        // }

      }

    }

    flmt = flmt(arma::find(flmt != id0));

  }

  arma::uvec f(m.n_cols, arma::fill::zeros); f(id) = 1;

  return gbp3q(p, fit.it, m, fit.k, f, fit.o, fit.ok);
}

// //' assert_matrix_last_col_dominant
// //' @description
// //'  auxilium of gbp3d_solver_dpp_filt
// //' @details
// //'  assert matrix last column dominant all other column: arma::all(m.col(i) <= m.col(m.n_cols - 1)), all i.
// //' @return ok_last_col_dominant
// //' @family gbp3d_solver_dpp_filt
// //' @export
// // [[Rcpp::export]]
// bool assert_matrix_last_col_dominant(const arma::mat& m) {
//
//   bool ok_last_col_dominant = true;
//
//   for (arma::uword i = 0; i < m.n_cols - 1; i++) {
//     // should make sure m such that l >= d >= h or vice versa
//     if (arma::all(m.col(i) <= m.col(m.n_cols - 1) + tol)) { ok_last_col_dominant = false; break; }
//   }
//
//   return ok_last_col_dominant;
// }

// //' gbp3d_solver_dpp_filt_fast
// //' @description
// //'  auxilium of gbp3d_solver_dpp_filt
// //' @details
// //'  gbp3d_solver_dpp_filt_fast use in context of gbp3d_solver_dpp_filt:
// //'   when determine all it can fit into some one bn, then can safely filt
// //'   out bn such that cannot hold all it at once determined via comparing
// //'   scale and volume
// //' @return flmt
// //' @family gbp3d_solver_dpp_filt
// //' @note
// //'  call in r won't see result since flmt is changed in place - designed to use in rcpp
// //' @export
// // [[Rcpp::export]]
void gbp3d_solver_dpp_filt_fast(
  const arma::mat& sldh, const double vldh, const arma::mat& sm, const arma::rowvec& vm, arma::uvec& flmt
) {

  // main
  arma::uvec flmtidk(flmt.size(), arma::fill::zeros);

  for (arma::uword i = 0; i < flmt.size(); i++) {
    if (arma::any(sldh.row(0) > sm(0, flmt(i))) ||
        arma::any(sldh.row(1) > sm(1, flmt(i))) ||
        arma::any(sldh.row(2) > sm(2, flmt(i))) ||
        vldh > vm(flmt(i))
    ) {
      flmtidk(i) = 1;
    }
  }

  flmt = flmt(arma::find(flmtidk != 1));

}

// //' gbp3d_solver_dpp_filt_slow
// //' @description
// //'  auxilium of gbp3d_solver_dpp_filt
// //' @details
// //' gbp3d_solver_dpp_filt_slow use in context of gbp3d_solver_dpp_filt:
// //'  when determine some it with volume v0 can fit into some one bn,
// //'  can safely filt out bn such that the bn volume is less than v0.
// //' @return flmt
// //' @family gbp3d_solver_dpp_filt
// //' @note
// //'  call in r won't see result since flmt is changed in place - designed to use in rcpp
// //' @export
// // [[Rcpp::export]]
void gbp3d_solver_dpp_filt_slow(
  const double vldh, const arma::rowvec& vm, arma::uvec& flmt
) {

  // main
  flmt = flmt(arma::find(vm(flmt) >= vldh));

}

// //' gbp3d_solver_dpp_filt_dcol
// //' @description
// //'  auxilium of gbp3d_solver_dpp_filt
// //' @details
// //'  gbp3d_solver_dpp_filt_dcol use in context of gbp3d_solver_dpp_filt:
// //'   when determine all it can fit into some one bn, and a particular bn bn0
// //'   cannot fit all it, then can safely filt out bn1 dominated by bn0 in the
// //'   sense that bn1 also cannot fit all it.
// //'  bn1 dominated by bn0 if all bn1(l, d, h) <= bn0(l, d, h), l >= d >= h.
// //'  must first determine all it can fit into some one bn, otherwise it is possible
// //'   that no bn can fit all it, so even it is correct that if bn0 cannot fit all,
// //'   then bn1 cannot neither, but bn1 is smaller than bn0, and could fit the same
// //'   set of it and thus provide a higher preference w.r.t utilization rate.
// //' @return flmt
// //' @family gbp3d_solver_dpp_filt
// //' @note
// //'  call in r won't see result since flmt is changed in place - designed to use in rcpp
// //' @export
// // [[Rcpp::export]]
void gbp3d_solver_dpp_filt_dcol(
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
