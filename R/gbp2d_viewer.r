#------------------------------------------------------------------------------#
#------------------------------- gbp2d::gbp2d.r -------------------------------#
#------------------------- author: gyang274@gmail.com -------------------------#
#------------------------------------------------------------------------------#

#--------+---------+---------+---------+---------+---------+---------+---------#
#234567890123456789012345678901234567890123456789012345678901234567890123456789#

#------------------------------------------------------------------------------#
#------------------------------------ main ------------------------------------#
#------------------------------------------------------------------------------#

#' gbp2d_viewer
#'
#' @description
#'
#'  gbp2d solution viewer
#'
#' @param sn
#'
#'  sn gbp2d object, solution from gbp2d_solver_dpp, see gbp2d.
#'
#' @param title
#'
#'  title <character>
#'
#' @param subtitle
#'
#'  subtitle <character>
#'
#' @family gbp2d_viewer
#' @rdname gbp2d_viewer
#' @export
gbp2d_viewer <- function(
  sn, title = NULL, subtitle = NULL
) {

  if (!(class(sn) == "Rcpp_gbp2d")) {
    stop("gbp2d_viewer: sn must be a gbp2d object with class Rcpp_gbp2d.")
  }

  kk <- which(sn$k == 1)

  it <- data.table(
    id = kk,
    x  = sn$it[1, kk],
    y  = 0, # stereo perception
    z  = sn$it[2, kk],
    l  = sn$it[3, kk],
    d  = 1, # stereo perception
    h  = sn$it[4, kk]
  )

  bn <- data.table(
    id = 1,
    l  = sn$bn[1],
    d  = 1,
    h  = sn$bn[2]
  )

  if (is.null(title)) {
    title <- "gbp2d_viewer"
  }

  if (is.null(subtitle)) {
    subtitle <- paste0("bin: (", bn[["l"]], ", ", bn[["h"]], "): # fit: ", length(kk), " of ", length(sn$k))
  }

  bpp_viewer_single(
    it, bn, title = title, subtitle = subtitle
  )

}


#' gbp2q_viewer
#'
#' @description
#'
#'  gbp2q solution viewer
#'
#' @param sn
#'
#'  sn gbp2q object, solution from gbp2d_solver_dpp_filt, see gbp2q.
#'
#' @param title
#'
#'  title <character>
#'
#' @param subtitle
#'
#'  subtitle <character>
#'
#' @family gbp2q_viewer
#' @rdname gbp2q_viewer
#' @export
gbp2q_viewer <- function(
  sn, title = NULL, subtitle = NULL
) {

  if (!(class(sn) == "Rcpp_gbp2q")) {
    stop("gbp2q_viewer: sn must be a gbp2q object with class Rcpp_gbp2q.")
  }

  kk <- which(sn$k == 1)

  it <- data.table(
    id = kk,
    x  = sn$it[1, kk],
    y  = 0,
    z  = sn$it[2, kk],
    l  = sn$it[3, kk],
    d  = 1,
    h  = sn$it[4, kk]
  )

  ff <- which(sn$f == 1)

  bn <- data.table(
    id = ff,
    l = sn$bn[1, ff],
    d = 1,
    h = sn$bn[2, ff]
  )

  if (is.null(title)) {
    title <- "gbp2q_viewer"
  }

  if (is.null(subtitle)) {
    subtitle <- paste0("bin: (", bn[["l"]], ", ", bn[["h"]], "): # fit: ", length(kk), " of ", length(sn$k))
  }

  bpp_viewer_single(
    it, bn, title = title, subtitle = subtitle
  )

}

#------------------------------------------------------------------------------#
