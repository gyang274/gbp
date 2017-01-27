#ifndef __BINPACK_BINPACK_M__
#define __BINPACK_BINPACK_M__


#include "gbp_u.h"
#include "gbp1d.h"
#include "gbp2d_it.h"
#include "gbp2d.h"
#include "gbp2d_ck.h"
#include "gbp3d_it.h"
#include "gbp3d.h"
#include "gbp3d_ck.h"
#include "gbp4d_it.h"
#include "gbp4d.h"
#include "gbp4d_ck.h"
#include "bpp.h"


// rcpp class and module -> r

// gbp1d.h
RCPP_EXPOSED_CLASS(gbp1d);

RCPP_MODULE(gbp1d_module) {

  using namespace Rcpp;

  class_<gbp1d>( "gbp1d" )

    .constructor<arma::vec, arma::uvec, arma::uword, arma::uvec, double, bool>()

    .field( "p" , &gbp1d::p  )
    .field( "w" , &gbp1d::w  )
    .field( "c" , &gbp1d::c  )
    .field( "k" , &gbp1d::k  )
    .field( "o" , &gbp1d::o  )
    .field( "ok", &gbp1d::ok )
    ;

  function(
    "gbp1d_solver_dpp", &gbp1d_solver_dpp, Rcpp::List::create(_["p"], _["w"], _["c"]),
    "gbp1d gbp1d_solver_dpp(const arma::vec& p, const arma::uvec& w, const arma::uword c)"
  );

}

// gbp2d_it.h
RCPP_EXPOSED_CLASS(Ktlist2d);

RCPP_MODULE(Ktlist2d_module) {

  using namespace Rcpp;

  class_<Ktlist2d>( "Ktlist2d" )

    .constructor<arma::uword, arma::mat, arma::field<arma::mat>, arma::vec>()

    .field( "n" , &Ktlist2d::n  )
    .field( "kt", &Ktlist2d::kt )
    .field( "xp", &Ktlist2d::xp )
    .field( "s" , &Ktlist2d::s  )

    ;

  function(
    "gbp2d_it_create_ktlist", &gbp2d_it_create_ktlist, Rcpp::List::create(_["bn"], _["it"], _["xp"], _["ktinit"], _["nlmt"]),
    "Ktlist2d gbp2d_it_create_ktlist(const arma::vec& bn, const arma::mat& it, const arma::mat& xp, const arma::vec& ktinit, const arma::uword nlmt)"
  );

}

// gbp2d.h
RCPP_EXPOSED_CLASS(gbp2d);

RCPP_MODULE(gbp2d_module) {

  using namespace Rcpp;

  class_<gbp2d>( "gbp2d" )

    .constructor<arma::vec, arma::mat, arma::vec, arma::uvec, double, bool>()

    .field( "p" , &gbp2d::p  )
    .field( "it", &gbp2d::it )
    .field( "bn", &gbp2d::bn )
    .field( "k" , &gbp2d::k  )
    .field( "o" , &gbp2d::o  )
    .field( "ok", &gbp2d::ok )

    ;

  function(
    "gbp2d_solver_dpp", &gbp2d_solver_dpp, Rcpp::List::create(_["p"], _["ld"], _["m"]),
    "gbp2d gbp2d_solver_dpp(const arma::vec& p, const arma::mat& ld, const arma::vec& m)"
  );

  function(
    "gbp2d_checkr", &gbp2d_checkr, Rcpp::List::create(_["sn"]), "bool gbp2d_checkr(gbp2d sn)"
  );

}

// gbp2d.h
RCPP_EXPOSED_CLASS(gbp2q);

RCPP_MODULE(gbp2q_module) {

  using namespace Rcpp;

  class_<gbp2q>( "gbp2q" )

    .constructor<arma::vec, arma::mat, arma::mat, arma::uvec, arma::uvec, double, bool>()

    .field( "p" , &gbp2q::p  )
    .field( "it", &gbp2q::it )
    .field( "bn", &gbp2q::bn )
    .field( "k" , &gbp2q::k  )
    .field( "f" , &gbp2q::f  )
    .field( "o" , &gbp2q::o  )
    .field( "ok", &gbp2q::ok )

    ;

  function(
    "gbp2d_solver_dpp_filt", &gbp2d_solver_dpp_filt, Rcpp::List::create(_["ld"], _["m"]),
    "gbp2q gbp2d_solver_dpp_filt(const arma::mat& ld, const arma::mat& m)"
  );

  function(
    "gbp2q_checkr", &gbp2q_checkr, Rcpp::List::create(_["sn"]), "bool gbp2q_checkr(gbp2q sn)"
  );

}

// gbp3d_it.h
RCPP_EXPOSED_CLASS(Ktlist3d);

RCPP_MODULE(Ktlist3d_module) {

  using namespace Rcpp;

  class_<Ktlist3d>( "Ktlist3d" )

    .constructor<arma::uword, arma::mat, arma::field<arma::mat>, arma::vec>()

    .field( "n" , &Ktlist3d::n  )
    .field( "kt", &Ktlist3d::kt )
    .field( "xp", &Ktlist3d::xp )
    .field( "s" , &Ktlist3d::s  )

    ;

  function(
    "gbp3d_it_create_ktlist", &gbp3d_it_create_ktlist, Rcpp::List::create(_["bn"], _["it"], _["xp"], _["ktinit"], _["nlmt"]),
    "Ktlist3d gbp3d_it_create_ktlist(const arma::vec& bn, const arma::mat& it, const arma::mat& xp, const arma::vec& ktinit, const arma::uword nlmt)"
  );

}

// gbp3d.h
RCPP_EXPOSED_CLASS(gbp3d);

RCPP_MODULE(gbp3d_module) {

  using namespace Rcpp;

  class_<gbp3d>( "gbp3d" )

    .constructor<arma::vec, arma::mat, arma::vec, arma::uvec, double, bool>()

    .field( "p" , &gbp3d::p  )
    .field( "it", &gbp3d::it )
    .field( "bn", &gbp3d::bn )
    .field( "k" , &gbp3d::k  )
    .field( "o" , &gbp3d::o  )
    .field( "ok", &gbp3d::ok )

    ;

  function(
    "gbp3d_solver_dpp", &gbp3d_solver_dpp, Rcpp::List::create(_["p"], _["ldh"], _["m"]),
    "gbp3d gbp3d_solver_dpp(const arma::vec& p, const arma::mat& ldh, const arma::vec& m)"
  );

  function(
    "gbp3d_checkr", &gbp3d_checkr, Rcpp::List::create(_["sn"]), "bool gbp3d_checkr(gbp3d sn)"
  );

}

// gbp3d.h
RCPP_EXPOSED_CLASS(gbp3q);

RCPP_MODULE(gbp3q_module) {

  using namespace Rcpp;

  class_<gbp3q>( "gbp3q" )

    .constructor<arma::vec, arma::mat, arma::mat, arma::uvec, arma::uvec, double, bool>()

    .field( "p" , &gbp3q::p  )
    .field( "it", &gbp3q::it )
    .field( "bn", &gbp3q::bn )
    .field( "k" , &gbp3q::k  )
    .field( "f" , &gbp3q::f  )
    .field( "o" , &gbp3q::o  )
    .field( "ok", &gbp3q::ok )

    ;

  function(
    "gbp3d_solver_dpp_filt", &gbp3d_solver_dpp_filt, Rcpp::List::create(_["ldh"], _["m"]),
    "gbp3q gbp3d_solver_dpp_filt(const arma::mat& ldh, const arma::mat& m)"
  );

  function(
    "gbp3q_checkr", &gbp3q_checkr, Rcpp::List::create(_["sn"]), "bool gbp3q_checkr(gbp3q sn)"
  );

}

// gbp4d_it.h
RCPP_EXPOSED_CLASS(Ktlist4d);

RCPP_MODULE(Ktlist4d_module) {

  using namespace Rcpp;

  class_<Ktlist4d>( "Ktlist4d" )

    .constructor<arma::uword, arma::mat, arma::field<arma::mat>, arma::vec>()

    .field( "n" , &Ktlist4d::n  )
    .field( "kt", &Ktlist4d::kt )
    .field( "xp", &Ktlist4d::xp )
    .field( "s" , &Ktlist4d::s  )

    ;

  function(
    "gbp4d_it_create_ktlist", &gbp4d_it_create_ktlist, Rcpp::List::create(_["bn"], _["it"], _["xp"], _["ktinit"], _["nlmt"]),
    "Ktlist4d gbp4d_it_create_ktlist(const arma::vec& bn, const arma::mat& it, const arma::mat& xp, const arma::vec& ktinit, const arma::uword nlmt)"
  );

}

// gbp4d.h
RCPP_EXPOSED_CLASS(gbp4d);

RCPP_MODULE(gbp4d_module) {

  using namespace Rcpp;

  class_<gbp4d>( "gbp4d" )

    .constructor<arma::vec, arma::mat, arma::vec, arma::uvec, double, bool>()

    .field( "p" , &gbp4d::p  )
    .field( "it", &gbp4d::it )
    .field( "bn", &gbp4d::bn )
    .field( "k" , &gbp4d::k  )
    .field( "o" , &gbp4d::o  )
    .field( "ok", &gbp4d::ok )

    ;

  function(
    "gbp4d_solver_dpp", &gbp4d_solver_dpp, Rcpp::List::create(_["p"], _["ldhw"], _["m"]),
    "gbp4d gbp4d_solver_dpp(const arma::vec& p, const arma::mat& ldhw, const arma::vec& m)"
  );

  function(
    "gbp4d_checkr", &gbp4d_checkr, Rcpp::List::create(_["sn"]), "bool gbp4d_checkr(gbp4d sn)"
  );

}

// gbp4d.h
RCPP_EXPOSED_CLASS(gbp4q);

RCPP_MODULE(gbp4q_module) {

  using namespace Rcpp;

  class_<gbp4q>( "gbp4q" )

    .constructor<arma::vec, arma::mat, arma::mat, arma::uvec, arma::uvec, double, bool>()

    .field( "p" , &gbp4q::p  )
    .field( "it", &gbp4q::it )
    .field( "bn", &gbp4q::bn )
    .field( "k" , &gbp4q::k  )
    .field( "f" , &gbp4q::f  )
    .field( "o" , &gbp4q::o  )
    .field( "ok", &gbp4q::ok )

    ;

  function(
    "gbp4d_solver_dpp_filt", &gbp4d_solver_dpp_filt, Rcpp::List::create(_["ldhw"], _["m"]),
    "gbp4q gbp4d_solver_dpp_filt(const arma::mat& ldhw, const arma::mat& m)"
  );

  function(
    "gbp4q_checkr", &gbp4q_checkr, Rcpp::List::create(_["sn"]), "bool gbp4q_checkr(gbp4q sn)"
  );

}

// bpp.h
RCPP_EXPOSED_CLASS(bppSgl);

RCPP_MODULE(bppSgl_module) {

  using namespace Rcpp;

  class_<bppSgl>( "bppSgl" )

    .constructor<arma::uvec, arma::mat, arma::mat, arma::uvec, arma::uvec, bool>()

    .field( "id", &bppSgl::id )
    .field( "it", &bppSgl::it )
    .field( "bn", &bppSgl::bn )
    .field( "k" , &bppSgl::k  )
    .field( "kb", &bppSgl::kb )
    .field( "ok", &bppSgl::ok )

    ;

  function(
    "bpp_solver_dpp", &bpp_solver_dpp, Rcpp::List::create(_["id"], _["ldhw"], _["m"]),
    "bppSgl bpp_solver_dpp(const arma::uvec& id, const arma::mat& ldhw, const arma::mat& m)"
  );

  function(
    "bpp_solver_sgl", &bpp_solver_sgl, Rcpp::List::create(_["ldhw"], _["m"]),
    "bppSgl bpp_solver_sgl(const arma::mat& ldhw, const arma::mat& m)"
  );

}


#endif // __BINPACK_BINPACK_M__
