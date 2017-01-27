#ifndef __BINPACK_BINPACK2D_XP__
#define __BINPACK_BINPACK2D_XP__


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// using namespace arma;


// gbp2d_xp.cpp

arma::mat gbp2d_xp_create_xp(const arma::vec& bn, const arma::mat& it);

void gbp2d_xp_update_xp(const arma::vec& bn, const arma::mat& it, const arma::vec& kt, arma::mat& xp);

void gbp2d_xp_purify_xp(arma::mat& xp);

void gbp2d_xp_update_xp_ikt(const arma::mat& it, const arma::vec& kt, arma::mat& xp);

void gbp2d_xp_update_rs_spg(const arma::mat& it, const arma::vec& kt, arma::mat& minBound, arma::mat& xpUpdate);

void gbp2d_xp_update_minbnd(const arma::vec& it, const arma::vec& kt, arma::mat& minBound, arma::mat& xpUpdate);

void gbp2d_xp_update_xp_spg(const arma::mat& it, const arma::vec& kt, arma::vec& maxBound, arma::mat& xpUpdate);

void gbp2d_xp_update_maxbnd(const arma::vec& it, const arma::vec& kt, arma::vec& maxBound, arma::mat& xpUpdate);

arma::uvec gbp2d_xp_it_qjt_kt(const arma::vec& it, const arma::vec& kt);

arma::uvec gbp2d_xp_it_pjt_kt(const arma::vec& it, const arma::vec& kt);


#endif // __BINPACK_BINPACK2D_XP__
