#ifndef __BINPACK_BINPACK1D__
#define __BINPACK_BINPACK1D__


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// using namespace arma;


// gbp1d
// @description
//  generalized bin packing problem in 1 dimension, a.k.a knapsack 0-1 problem.
// @details
//  gbp1d init a profit vector p, a weight vector w, and a weight constraint c,
//  gbp1d solver would solve
//
//    maximize   sum_{j=1}^{n} p_{j} k_{j}
//
//    subject to sum_{j=1}^{n} w_{j} k_{j} leq c
//               k_{j} in {0, 1}, j = 1 , .. , n
//
//  and instantiate a gbp1d object with a selectin vector k and an objective o.
// @family gbp1d
// @rdname gbp1d
class gbp1d {

public:

  arma::vec p; arma::uvec w; arma::uword c; arma::uvec k; double o; bool ok;

  // .constructor
  gbp1d(arma::vec p, arma::uvec w, arma::uword c, arma::uvec k, double o, bool ok) : p(p), w(w), c(c), k(k), o(o), ok(ok) {}

};


// solve gbp1d via dynamic programming simple - adagio::knapsnak()
gbp1d gbp1d_solver_dpp(const arma::vec& p, const arma::uvec& w, const arma::uword c);

// solve gbp1d via dynamic programming w. minimized search states - minknap.c
// gbp1d gbp3d_solver_min(const arma::vec& p, const arma::uvec& w, const arma::uword c);


#endif // __BINPACK_BINPACK1D__
