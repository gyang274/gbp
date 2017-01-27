#' gbp1d
#' @aliases
#'  gbp1d Rcpp_gbp1d Rcpp_gbp1d-class
#' @description
#'  generalized bin packing problem in 1 dimension, a.k.a knapsack 0-1 problem.
#' @details
#'  gbp1d init a profit vector p, a weight vector w, and a weight constraint c,
#'  gbp1d solver would solve
#'
#'     maximize sum_{j=1}^{n} p_{j} x_{j}
#'
#'   subject to sum_{j=1}^{n} w_{j} x_{j} leq c
#'              x_{j} in {0, 1}, j = 1, ...., n
#'
#'  and instantiate a gbp1d object with a selectin vector x and an objective z.
#'
#'  gbp1d is implemented as rcpp class, an instantiate can be solved by calling
#'   gbp1d_solver_dpp(p, w, c) and gbp1d_solver_min(p, w, c)
#' @family gbp1d
#' @rdname gbp1d
#' @docType class
"gbp1d"

#' gbp1d_solver_dpp
#' @description
#'  solve gbp1d via dynamic programming simple - adagio::knapsnak()
#' @details
#'  a dynamic programming solver on gbp1d instantiate - knapsack 0-1 problem, see gbp1d.
#'
#'  gbp1d init a profit vector p, a weight vector w, and a weight constraint c,
#'  gbp1d solver would solve
#'
#'     maximize sum_{j=1}^{n} p_{j} x_{j}
#'
#'   subject to sum_{j=1}^{n} w_{j} x_{j} leq c
#'              x_{j} in {0, 1}, j = 1, ...., n
#'
#'  and instantiate a gbp1d object with a selectin vector x and an objective z.
#'
#'  gbp1d is implemented as rcpp class, an instantiate can be solved by calling
#'   gbp1d_solver_dpp(p, w, c) and gbp1d_solver_min(p, w, c)
#'
#' @param p
#'  p profit <vector>::<numeric>
#' @param w
#'  w weight <vector>::<integer>
#' @param c
#'  c constraint on weight <integer>
#' @return gbp1d
#'  a gbp1d instantiate with p profit, w weight, c constraint on weight,
#'   k selection, o objective, and ok an indicator of all fit or not.
#' @family gbp1d
#' @rdname gbp1d_solver_dpp
"gbp1d_solver_dpp"
