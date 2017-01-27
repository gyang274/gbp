#ifndef __BINPACK_BINPACK_BPP__
#define __BINPACK_BINPACK_BPP__


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// using namespace arma;


// bppSgl
// @description
//  bpp solution of a single one order or multiple order
// @details
//  packing it into multiple bn w.r.t bn size and weight limit while select bn as small as possible
// @param id <vector>
//  id order id <integer> list - should sorted or at least grouped w.r.t order id
// @param it
//  it position and scale <matrix>
//  - x, y, z, w it position and w hold in the bin <numeric>
//  - l, d, h, w it scale along x, y, z and also w <numeric>
// @param bn
//  bn scale <vector>
//  - l, d, h, w bn scale along x, y, z and also w <numeric>
// @param k
//  k ticket id indicator 0 (if cannot fit into any bin), 1, 2, 3, 4, ... <vector>
// @param kb
//  kb ticket bn id indicator - which bn to use for packing each ticket <vector>
// @param ok
//  ok a quick indicator of any it can not fit into any bn? <bool>
class bppSgl {

public:

  arma::uvec id; arma::mat it; arma::mat bn; arma::uvec k; arma::uvec kb; bool ok;

  // .constructor
  bppSgl(
   arma::uvec id, arma::mat it, arma::mat bn, arma::uvec k, arma::uvec kb, bool ok
  ) : id(id), it(it), bn(bn), k(k), kb(kb), ok(ok) {
  }

};

bppSgl bpp_solver_dpp(
  const arma::uvec& id, const arma::mat& ldhw, const arma::mat& m
);

bppSgl bpp_solver_sgl(
  const arma::mat& ldhw, const arma::mat& m
);

bool bpp_solver_sgl_screen(
  const arma::mat& ldhw, const arma::mat& m, arma::uvec& itlmt
);

#endif // __BINPACK_BINPACK4D__
