#------------------------------------------------------------------------------#
#--------------------------------- gbp::gbp.r ---------------------------------#
#------------------------- author: gyang274@gmail.com -------------------------#
#------------------------------------------------------------------------------#

#--------+---------+---------+---------+---------+---------+---------+---------#
#234567890123456789012345678901234567890123456789012345678901234567890123456789#

#------------------------------------------------------------------------------#
#------------------------------------ main ------------------------------------#
#------------------------------------------------------------------------------#

#' bpp_solver
#'
#' @description
#'
#'  bpp single or multiple order packing solver
#'
#' @details
#'
#'  bpp solver is designed to solve packing in warehouse
#'
#'  bpp solver digest input it as a list of order (oid) and each row contain one
#'   sku (sku) in an order with length (l), depth (d), height (h) and weight (w)
#'   and aims to pack it list into one or more bin from a given list of bin that
#'   bin length (l), depth (d), height (h), and a single weight limit (wlmt).
#'
#'  bn list must be sorted by volume so that the smaller the eariler and preferred,
#'   and each bn must be sorted so that l >= d >= h
#'
#'  bpp solver would call bpp_solver_dpp_wrapper and aims to find a packing schema
#'   such that: use as small number of bin as possible, and use small bin whenever
#'   possible, w.r.t the 3d none overlap constraint and weight limit constraint.
#'
#' @param it
#'
#'  it item <data.table>
#'
#'  - oid: order id <integer>
#'
#'  - sku: stock keeping unit as it id <character>
#'
#'  - l: it length which scale will be placed along x-coordinate <numeric>
#'
#'  - d: it depth  which scale will be placed along y-coordinate <numeric>
#'
#'  - h: it height which scale will be placed along z-coordinate <numeric>
#'
#'  - w: it weight optional which scale will be used restriction <integer>
#'
#' @param bn
#'
#'  bn bins <data.table>
#'
#'  - id: bn id <character>
#'
#'  - l: bn length limit along x-coordinate <numeric>
#'
#'  - d: bn depth  limit along y-coordinate <numeric>
#'
#'  - h: bn height limit along z-coordinate <numeric>
#'
#'  - w: bn weight limit along w - a separate single dimension <numeric>
#'
#'  - l, d, h will be sorted to have l >= d >= h within solver
#'
#' @return sn
#'
#'  sn bpp_solution <list>
#'
#'  - it item <data.table>
#'
#'    - oid: order id <integer>
#'
#'    - sku: stock keeping unit as it id <character>
#'
#'    - tid: ticket id - an unique id within oid <integer>
#'
#'    - otid: order id x ticket id - an unique indentifier indicate it with same tid can be packed into one bin <character>
#'
#'    - bid: bn id <integer>
#'
#'    - x, y, z it position in bid bin <numeric>
#'
#'    - l, d, h it scale along x, y, z <numeric>
#'
#'    - w it weight <numeric>
#'
#'  - bn bins <data.table>
#'
#'    - id bn id <character>
#'
#'    - l bn length limit along x-coordinate <numeric>
#'
#'    - d bn depth  limit along y-coordinate <numeric>
#'
#'    - h bn height limit along z-coordinate <numeric>
#'
#'    - w bn weight limit along w - a separate single dimension <numeric>
#'
#' @note
#'  bpp_solver is an r-level wrapper over c-level bpp_solver_dpp_wrapper,
#'   add otid as an unique indentifier.
#' @family bpp_solver
#' @rdname bpp_solver
#' @export
bpp_solver <- function(it, bn) {

  #- init
  lkit <- copy(it) %>% `class<-`(c("data.table", "data.frame"))

  lkit <- lkit %>% data.table::setorderv(c("oid", "sku"))

  lkit <- lkit[ , list(
    oid = get("oid"),
    sku = get("sku"),
    l   = get("l"),
    d   = get("d"),
    h   = get("h"),
    w   = w
  )]

  #- init
  lkbn <- copy(bn) %>% `class<-`(c("data.table", "data.frame"))

  lkbn <- lkbn[ , list(
    id = get("id"),
    l  = pmax(get("l"), get("d"), get("h")),
    d  = apply(cbind(get("l"), get("d"), get("h")), 1L, stats::median),
    h  = pmin(get("l"), get("d"), get("h")),
    w  = get("w")
  )]

  #- main
  lkit <- gbp::bpp_solver_dpp_wrapper(lkit, lkbn) %>% `class<-`(c("data.table", "data.frame"))

  lkit <- copy(lkit)

  lkit[ , `:=`(otid = paste0(get("oid"), "X", get("tid")))]

  data.table::setcolorder(lkit, c("oid", "tid", "otid", "bid", "sku", "x", "y", "z", "l", "d", "h", "w"))

  return(list(it = lkit, bn = lkbn))

}
#------------------------------------------------------------------------------#

