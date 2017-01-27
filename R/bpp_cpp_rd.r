#' bppSgl
#' @aliases
#'  bppSgl Rcpp_bppSgl Rcpp_bppSgl-class
#' @description
#'  bpp solution of a single one order or multiple order
#' @details
#'  packing it into multiple bn w.r.t bn size and weight limit while select bn as small as possible
#'
#'  a bppSgl class instance has 6 fields:
#'
#'   - id: order id <integer>
#'
#'      list - should sorted or at least grouped w.r.t order id
#'
#'   - it: it position and scale <matrix>
#'
#'     - x, y, z, w it position and w in the bin <numeric> (w hold in bn when fit it in bn)
#'
#'     - l, d, h, w it scale along x, y, z and w <numeric> (w of it itself)
#'
#'   - bn: bn scale <vector>
#'
#'     - l, d, h, w bn scale along x, y, z and w <numeric>
#'
#'   - k: ticket id indicator 0 (if cannot fit into any bin), 1, 2, 3, 4, ... <vector>
#'
#'   - kb: ticket bn id indicator - which bn to use for packing each ticket <vector>
#'
#'   - ok: ok a quick indicator of any it can not fit into any bn? <bool>
#'
#' @family bppSgl
#' @rdname bppSgl
#' @docType class
#' @export
"bppSgl"


#' bpp_solver_dpp
#'
#' @description
#'
#'  main solver of e-commerce warehouse packing algorithm
#'
#' @details
#'
#'  bpp init a list of order on sku in data.frame it - oid, sku, l, d, h, w:
#'   order id oid, stock keeping unit sku, length l, depth d, height h and weight w,
#'
#'  and also a list of available bn in data.frame bn - id, l, d, h, w:
#'   bn id, length l, depth d, height h, and weight limit w, sorted by peference often smaller prefered,
#'
#'  and a single weight limit wlmt applied on all bin.
#'
#'  bpp solver would solve
#'
#'   select least number of bn for packing each order w.r.t bn size and weight limit
#'    and make sure the bn selected are as small as possible.
#'
#' @param id <vector>
#'
#'  id order id <integer> vector - should sorted or at least grouped w.r.t order id
#'
#' @param ldhw <matrix>
#'
#'  it order list
#'
#'  - l, d, h, w it scale along x, y, z and w <numeric>
#'
#'  it columns should corresponding to id
#'
#' @param m <matrix>
#'
#'  m a bin list
#'
#'  - l, d, h, w bn scale along x, y, z and w <numeric>
#'
#'  m should sorted w.r.t preference
#'
#' @return bppSgl
#' @family bpp_solver_dpp
#' @export
"bpp_solver_dpp"


#' bpp_solver_sgl
#' @description
#'  subroutine of bpp_solver_dpp
#' @details
#'  fit a single order into bn list, call gbp4d_solver_dpp_filt() as main solver.
#' @inheritParams bpp_solver_dpp
#' @return bppSgl
#' @family bpp_solver_dpp
#' @export
"bpp_solver_sgl"
