#------------------------------------------------------------------------------#
#------------------------------- gbp4d::gbp4d.r -------------------------------#
#------------------------- author: gyang274@gmail.com -------------------------#
#------------------------------------------------------------------------------#

#--------+---------+---------+---------+---------+---------+---------+---------#
#234567890123456789012345678901234567890123456789012345678901234567890123456789#

#------------------------------------------------------------------------------#
#------------------------------------ main ------------------------------------#
#------------------------------------------------------------------------------#

#' gbp4d_viewer
#'
#' @description
#'
#'  gbp4d solution viewer
#'
#' @param sn
#'
#'  sn gbp4d object, solution from gbp4d_solver_dpp, see gbp4d.
#'
#' @param title
#'
#'  title <character>
#'
#' @param subtitle
#'
#'  subtitle <character>
#'
#' @family gbp4d_viewer
#' @rdname gbp4d_viewer
#' @export
gbp4d_viewer <- function(
  sn, title = NULL, subtitle = NULL
) {

  if (!(class(sn) == "Rcpp_gbp4d")) {
    stop("gbp4d_viewer: sn must be a gbp4d object with class Rcpp_gbp4d.")
  }

  kk <- which(sn$k == 1)

  it <- data.table(
    id = kk,
    x  = sn$it[1, kk],
    y  = sn$it[2, kk],
    z  = sn$it[3, kk],
    # w  = sn$it[4, kk], # w hold in bin when fit it
    l  = sn$it[5, kk],
    d  = sn$it[6, kk],
    h  = sn$it[7, kk],
    w  = sn$it[8, kk] # w of it
  )

  bn <- data.table(
    id = 1,
    l  = sn$bn[1],
    d  = sn$bn[2],
    h  = sn$bn[3],
    w  = sn$bn[4]
  )

  if (is.null(title)) {
    title <- "gbp4d_viewer"
  }

  if (is.null(subtitle)) {
    subtitle <- paste0("bin: (", bn[["l"]], ", ", bn[["d"]], ", ", bn[["h"]], "): # fit: ", length(kk), " of ", length(sn$k), " - total weight:", sum(it[["w"]]))
  }

  bpp_viewer_single(
    it, bn, title = title, subtitle = subtitle
  )

}


#' gbp4q_viewer
#'
#' @description
#'
#'  gbp4q solution viewer
#'
#' @param sn
#'
#'  sn gbp4q object, solution from gbp4d_solver_dpp_filt, see gbp4q.
#'
#' @param title
#'
#'  title <character>
#'
#' @param subtitle
#'
#'  subtitle <character>
#'
#' @family gbp4q_viewer
#' @rdname gbp4q_viewer
#' @export
gbp4q_viewer <- function(
  sn, title = NULL, subtitle = NULL
) {

  if (!(class(sn) == "Rcpp_gbp4q")) {
    stop("gbp4q_viewer: sn must be a gbp4q object with class Rcpp_gbp4q.")
  }

  kk <- which(sn$k == 1)

  it <- data.table(
    id = kk,
    x  = sn$it[1, kk],
    y  = sn$it[2, kk],
    z  = sn$it[3, kk],
    # w  = sn$it[4, kk], # w hold in bin when fit it
    l  = sn$it[5, kk],
    d  = sn$it[6, kk],
    h  = sn$it[7, kk],
    w  = sn$it[8, kk] # w of it
  )

  ff <- which(sn$f == 1)

  bn <- data.table(
    id = ff,
    l  = sn$bn[1, ff],
    d  = sn$bn[2, ff],
    h  = sn$bn[3, ff],
    w  = sn$bn[4, ff]
  )

  if (is.null(title)) {
    title <- "gbp4q_viewer"
  }

  if (is.null(subtitle)) {
    subtitle <- paste0("bin: (", bn[["l"]], ", ", bn[["d"]], ", ", bn[["h"]], "): # fit: ", length(kk), " of ", length(sn$k), " - total weight:", sum(it[["w"]]))
  }

  bpp_viewer_single(
    it, bn, title = title, subtitle = subtitle
  )

}

#------------------------------------------------------------------------------#
