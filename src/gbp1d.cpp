#include "gbp1d.h"


// solve gbp1d via dynamic programming simple - adagio::knapsnak()
// gbp1d gbp1d_solver_dpp(const arma::vec& p, const arma::uvec& w, const arma::uword c);
gbp1d gbp1d_solver_dpp(const arma::vec& p, const arma::uvec& w, const arma::uword c) {

  // init
  arma::uword n = w.size();

  // init
  arma::uvec k = arma::zeros<arma::uvec>(n); double o = 0.00; bool ok = false;

  // init
  arma::mat q(c + 1, n, arma::fill::zeros);
  arma::vec g = arma::zeros<arma::vec>(c + 1);

  // forward
  for (arma::uword i = 0; i < n; i++) {
    q.col(i) = g;
    for (arma::uword j = w(i); j < c + 1; j++) {
      // must go reverse otherwise double count
      g(c - j + w(i)) = std::max(g(c - j + w(i)), g(c - j) + p(i));
    }
  }
  double qmax = g(c);

  // backward
  double qptr = qmax;
  arma::uword j = c;
  for (arma::uword i = 0; i < n; i++) {
    if ( q(j, n - 1 - i) < qptr ) {
      k(n - 1 - i) = 1;
      j = j - w(n - 1 - i);
      qptr = q(j, n - 1 - i);
    }
  }

  o = sum(p(arma::find(k == 1))); ok = arma::all(k == 1);

  return gbp1d(p, w, c, k, o, ok);
}


// solve gbp1d via dynamic programming w. search on minimized num core state - minknap.c
// gbp1d gbp1d_solver_min(const arma::vec& p, const arma::uvec& w, const arma::uword c);
