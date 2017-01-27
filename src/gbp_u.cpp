#include "gbp_u.h"


// //' create_max_mode_tbl
// //' @description
// //'  create a max mode tbl from a matrix
// //' @details
// //'  create_max_mode_tbl init a matrix m, a vector mrow ulmt, a vector mcol vlmt,
// //'  a max value constraint mval, a max level constraint mlvl, and will create a
// //'  max mode tbl via
// //'
// //'  init a value count mapping vc
// //'
// //'   all v = m(i, j) if ulmt(i) == 0, vlmt(j) == 0 and v <= mval then vc[v]++
// //'
// //'  and return mlvl value v with increasing count c
// //'
// //'   v_{0} = max{v: vc[v] > 0}; mv(0, _) = (v_{0}, vc[v_{0}]);
// //'
// //'   v_{i} = max{v: vc[v_{i}] > vc[v_{i-1}]}; mv(i, _) = (v_{i}, vc[v_{i}]);
// //' @family gbp_u
// //' @rdname create_max_mode_tbl
// // [[Rcpp::export]]
arma::umat create_max_mode_tbl(
    arma::umat& m, arma::uvec& ulmt, arma::uvec& vlmt,
    const int mval, const int mlvl
) {

  std::map<int, int> vc; //  value count map

  for (arma::uword i = 0; i < m.n_rows; i++) {
    if (ulmt(i) == 0) {
      for (arma::uword j = 0; j < m.n_cols; j++) {
        if (vlmt(j) == 0) vc[m(i, j)]++;
      }
    }
  }

  arma::umat mv(2, mlvl, arma::fill::zeros);

  int i = 0; int k = 0;
  for (std::map<int, int>::reverse_iterator it = vc.rbegin(); it != vc.rend(); it++) {
    if (it->first <= mval && it->second > k) {
      mv(0, i) = it->first; mv(1, i) = it->second; i++; k = it->second;
    }
    if (i == mlvl) break;
  }

  return mv;
}


// //' unique_rows
// //' @description
// //'  unique rows in a matrix
// //' @family gbp_u_unique
// //' @rdname unique_rows
// // [[Rcpp::export]]
arma::mat unique_rows(const arma::mat& m) {

  arma::uvec ulmt = arma::zeros<arma::uvec>(m.n_rows);

  for (arma::uword i = 0; i < m.n_rows; i++) {
    for (arma::uword j = i + 1; j < m.n_rows; j++) {
      if (approx_equal_cpp(m.row(i), m.row(j))) { ulmt(j) = 1; break; }
    }
  }

  return m.rows(find(ulmt == 0));
}

// //' unique_cols
// //' @description
// //'  unique cols in a matrix
// //' @family gbp_u_unique
// //' @rdname unique_cols
// // [[Rcpp::export]]
arma::mat unique_cols(const arma::mat& m) {

  arma::uvec vlmt = arma::zeros<arma::uvec>(m.n_cols);

  for (arma::uword i = 0; i < m.n_cols; i++) {
    for (arma::uword j = i + 1; j < m.n_cols; j++) {
      if (approx_equal_cpp(m.col(i), m.col(j))) { vlmt(j) = 1; break; }
    }
  }

  return m.cols(find(vlmt == 0));
}


// //' sort_index_via_cols_internal
// //' @description
// //'  sort a matrix via multiple cols
// //' @details
// //'  sort a matrix via specific cols
// //'   first sort rows in ulmt via col vlmt(0),
// //'   break ties rows in ulmt via col vlmt(1),
// //'   break ties rows in ulmt via col vlmt(2),
// //'   and so on and so on ..
// //'  internal use for recursive call
// //' @param m
// //'   matrix
// //' @param ulmt
// //'   rows indices in cpp style start from 0 - sort within ulmt rows
// //' @param vlmt
// //'   cols indices in cpp style start from 0 - sort ulmt rows with vlmt cols
// //' @return id
// //'   m matrix ulmt rows sorted indices in cpp style start from 0 sorted by vlmt cols
// //' @family gbp_u_sort
// //' @rdname sort_index_via_cols_internal
// // [[Rcpp::export]]
arma::uvec sort_index_via_cols_internal(const arma::mat& m, const arma::uvec& ulmt, const arma::uvec& vlmt) {

  if ( m.n_cols == 0 ) return arma::zeros<arma::uvec>(0);

  if ( vlmt.size() == 0 ) {
    if ( ulmt.size() == 0 ) {
      return arma::zeros<arma::uvec>(0);
    } else {
      return arma::linspace<arma::uvec>(0, ulmt.size() - 1, ulmt.size());
    }
  }

  arma::uvec mlmt = vlmt(0) * m.n_rows + ulmt;

  arma::uvec id = arma::stable_sort_index(m.elem(mlmt), "ascend");

  if ( vlmt.size() == 1 ) return id;

  if ( vlmt.size() > 1 ) {

    arma::uvec g(ulmt.size()); double v = arma::datum::nan;

    for (arma::uword i = 0, j = 0; i < ulmt.size(); i++) {
      if (v == m(ulmt(id(i)), vlmt(0))) {
        g(i) = j;
      } else {
        j = i; g(i) = j; v = m(ulmt(id(i)), vlmt(0));
      }
    }

    for (arma::uword i = 0; i < ulmt.size(); i++) {
      arma::uvec g0 = find(g == i);
      if (g0.size() > 1) {
        arma::uvec id0 = id(g0);
        arma::uvec id0id = sort_index_via_cols_internal(m, ulmt(id0), vlmt.subvec(1, vlmt.size() - 1));
        id(g0) = id0(id0id);
      }
    }

  }

  return id;
}

// //' sort_index_via_cols
// //' @description
// //'  sort a matrix via multiple cols
// //' @details
// //'  sort a matrix via specific cols
// //'   first sort via col vlmt(0),
// //'   break ties via col vlmt(1),
// //'   break ties via col vlmt(2),
// //'   and so on and so on ..
// //' @param m
// //'   matrix
// //' @param vlmt
// //'   cols indices in cpp style start from 0
// //' @return id
// //'   m matrix rows sorted indices in cpp style start from 0 sorted by vlmt cols
// //' @family gbp_u_sort
// //' @rdname sort_index_via_cols
// // [[Rcpp::export]]
arma::uvec sort_index_via_cols(const arma::mat& m, const arma::uvec& vlmt) {

  if ( m.n_cols == 0 || m.n_rows == 0 ) return arma::zeros<arma::uvec>(0);

  if ( vlmt.size() == 0 ) return arma::linspace<arma::uvec>(0, m.n_rows - 1, m.n_rows);

  arma::uvec id = sort_index_via_cols_internal(m, arma::linspace<arma::uvec>(0, m.n_rows - 1, m.n_rows), vlmt);

  return id;
}

// //' sort_via_cols
// //' @description
// //'  sort a matrix via multiple cols
// //' @details
// //'  sort a matrix via specific cols
// //'   first sort via col vlmt(0),
// //'   break ties via col vlmt(1),
// //'   break ties via col vlmt(2),
// //'   and so on and so on ..
// //' @inheritParams sort_index_via_cols
// //' @return ms
// //'   m matrix sorted by vlmt cols
// //' @family gbp_u_sort
// //' @rdname sort_via_cols
// // [[Rcpp::export]]
arma::mat sort_via_cols(const arma::mat& m, const arma::uvec& vlmt) {

  if ( m.n_cols == 0 || m.n_rows == 0 || vlmt.size() == 0 ) return m;

  arma::uvec id = sort_index_via_cols_internal(m, arma::linspace<arma::uvec>(0, m.n_rows - 1, m.n_rows), vlmt);

  return m.rows(id);
}


// //' sort_index_via_rows_internal
// //' @description
// //'  sort a matrix via multiple rows
// //' @details
// //'  sort a matrix via specific rows
// //'   first sort cols in vlmt via row ulmt(0),
// //'   break ties cols in vlmt via row ulmt(1),
// //'   break ties cols in vlmt via row ulmt(2),
// //'   and so on and so on ..
// //'  internal use for recursive call
// //' @param m
// //'   matrix
// //' @param ulmt
// //'   rows indices in cpp style start from 0 - sort vlmt cols with ulmt rows
// //' @param vlmt
// //'   cols indices in cpp style start from 0 - sort within vlmt cols
// //' @return id
// //'   m matrix vlmt cols sorted indices in cpp style start from 0 sorted by ulmt rows
// //' @family gbp_u_sort
// //' @rdname sort_index_via_rows_internal
// // [[Rcpp::export]]
arma::uvec sort_index_via_rows_internal(const arma::mat& m, const arma::uvec& ulmt, const arma::uvec& vlmt) {

  if ( m.n_rows == 0 ) return arma::zeros<arma::uvec>(0);

  if ( ulmt.size() == 0 ) {
    if ( vlmt.size() == 0 ) {
      return arma::zeros<arma::uvec>(0);
    } else {
      return arma::linspace<arma::uvec>(0, vlmt.size() - 1, vlmt.size());
    }
  }

  arma::uvec mlmt = vlmt * m.n_rows + ulmt(0);

  arma::uvec id = arma::stable_sort_index(m.elem(mlmt), "ascend");

  if ( ulmt.size() == 1 ) return id;

  if ( ulmt.size() > 1 ) {

    arma::uvec g(vlmt.size()); double v = arma::datum::nan;

    for (arma::uword i = 0, j = 0; i < vlmt.size(); i++) {
      if (v == m(ulmt(0), vlmt(id(i)))) {
        g(i) = j;
      } else {
        j = i; g(i) = j; v = m(ulmt(0), vlmt(id(i)));
      }
    }

    for (arma::uword i = 0; i < vlmt.size(); i++) {
      arma::uvec g0 = find(g == i);
      if (g0.size() > 1) {
        arma::uvec id0 = id(g0);
        arma::uvec id0id = sort_index_via_rows_internal(m, ulmt.subvec(1, ulmt.size() - 1), vlmt(id0));
        id(g0) = id0(id0id);
      }
    }

  }

  return id;
}

// //' sort_index_via_rows
// //' @description
// //'  sort a matrix via multiple rows
// //' @details
// //'  sort a matrix via specific rows
// //'   first sort via row ulmt(0),
// //'   break ties via row ulmt(1),
// //'   break ties via row ulmt(2),
// //'   and so on and so on ..
// //' @param m
// //'   matrix
// //' @param ulmt
// //'   rows indices in cpp style start from 0
// //' @return id
// //'   m matrix cols sorted indices in cpp style start from 0 sorted by ulmt rows
// //' @family gbp_u_sort
// //' @rdname sort_index_via_rows
// // [[Rcpp::export]]
arma::uvec sort_index_via_rows(const arma::mat& m, const arma::uvec& ulmt) {

  if ( m.n_rows == 0 || m.n_cols == 0 ) return arma::zeros<arma::uvec>(0);

  if ( ulmt.size() == 0 ) return arma::linspace<arma::uvec>(0, m.n_cols - 1, m.n_cols);

  arma::uvec id = sort_index_via_rows_internal(m, ulmt, arma::linspace<arma::uvec>(0, m.n_cols - 1, m.n_cols));

  return id;
}

// //' sort_via_rows
// //' @description
// //'  sort a matrix via multiple rows
// //' @details
// //'  sort a matrix via specific rows
// //'   first sort via row ulmt(0),
// //'   break ties via row ulmt(1),
// //'   break ties via row ulmt(2),
// //'   and so on and so on ..
// //' @inheritParams sort_index_via_rows
// //' @return ms
// //'   m matrix sorted by ulmt rows
// //' @family gbp_u_sort
// //' @rdname sort_via_rows
// // [[Rcpp::export]]
arma::mat sort_via_rows(const arma::mat& m, const arma::uvec& ulmt) {

  if ( m.n_rows == 0 || m.n_cols == 0 || ulmt.size() == 0 ) return m;

  arma::uvec id = sort_index_via_rows_internal(m, ulmt, arma::linspace<arma::uvec>(0, m.n_cols - 1, m.n_cols));

  return m.cols(id);
}

