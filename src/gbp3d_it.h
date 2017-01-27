#ifndef __BINPACK_BINPACK3D_IT__
#define __BINPACK_BINPACK3D_IT__


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// using namespace arma;


// gbp3d_it.cpp

// Ktlist3d
// @description
//  Ktlist3d hold multiple kt for recursive fit
// @detail
//  Ktlist3d hold multiple kt via consider all possible fit onto different xp and different rotation and nlimit
class Ktlist3d {

  public:

    arma::uword n; arma::mat kt; arma::field<arma::mat> xp; arma::vec s;

    // .constructor
    Ktlist3d(
      arma::uword n, arma::mat kt, arma::field<arma::mat> xp, arma::vec s
    ) : n(n), kt(kt), xp(xp), s(s) {
    }

    // // .constructor overload default n = 0 case - cannot export to r
    // Ktlist3d(
    // ) {
    //   n = 0; kt = arma::zeros<arma::mat>(6, 0); xp = arma::field<arma::mat>(0); s = arma::zeros<arma::vec>(0);
    // }

};

Ktlist3d gbp3d_it_create_ktlist(
  const arma::vec& bn, const arma::mat& it, const arma::mat& xp, const arma::vec& ktinit, const arma::uword nlmt
);

void gbp3d_it_purify_ktlist(arma::mat& kt, arma::field<arma::mat>& xplist, arma::vec& s, const arma::uword nlmt);

double gbp3d_it_scorer_ktlist(const arma::mat& xpWithKt);

arma::mat gbp3d_it_create_ktldht(const arma::vec& ktinit);


#endif // __BINPACK_BINPACK3D_IT__
