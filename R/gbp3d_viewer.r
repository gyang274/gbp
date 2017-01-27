#------------------------------------------------------------------------------#
#------------------------------- gbp3d::gbp3d.r -------------------------------#
#------------------------- author: gyang274@gmail.com -------------------------#
#------------------------------------------------------------------------------#

#--------+---------+---------+---------+---------+---------+---------+---------#
#234567890123456789012345678901234567890123456789012345678901234567890123456789#

#------------------------------------------------------------------------------#
#------------------------------------ main ------------------------------------#
#------------------------------------------------------------------------------#

#' gbp3d_viewer
#'
#' @description
#'
#'  gbp3d solution viewer
#'
#' @param sn
#'
#'  sn gbp3d object, solution from gbp3d_solver_dpp, see gbp3d.
#'
#' @param title
#'
#'  title <character>
#'
#' @param subtitle
#'
#'  subtitle <character>
#'
#' @family gbp3d_viewer
#' @rdname gbp3d_viewer
#' @export
gbp3d_viewer <- function(
  sn, title = NULL, subtitle = NULL
) {

  if (!(class(sn) == "Rcpp_gbp3d")) {
    stop("gbp3d_viewer: sn must be a gbp3d object with class Rcpp_gbp3d.")
  }

  kk <- which(sn$k == 1)

  it <- data.table(
    id = kk,
    x  = sn$it[1, kk],
    y  = sn$it[2, kk],
    z  = sn$it[3, kk],
    l  = sn$it[4, kk],
    d  = sn$it[5, kk],
    h  = sn$it[6, kk]
  )

  bn <- data.table(
    id = 1,
    l  = sn$bn[1],
    d  = sn$bn[2],
    h  = sn$bn[3]
  )

  if (is.null(title)) {
    title <- "gbp3d_viewer"
  }

  if (is.null(subtitle)) {
    subtitle <- paste0("bin: (", bn[["l"]], ", ", bn[["d"]], ", ", bn[["h"]], "): # fit: ", length(kk), " of ", length(sn$k))
  }

  bpp_viewer_single(
    it, bn, title = title, subtitle = subtitle
  )

}


#' gbp3q_viewer
#'
#' @description
#'
#'  gbp3q solution viewer
#'
#' @param sn
#'
#'  sn gbp3q object, solution from gbp3d_solver_dpp_filt, see gbp3q.
#'
#' @param title
#'
#'  title <character>
#'
#' @param subtitle
#'
#'  subtitle <character>
#'
#' @family gbp3q_viewer
#' @rdname gbp3q_viewer
#' @export
gbp3q_viewer <- function(
  sn, title = NULL, subtitle = NULL
) {

  if (!(class(sn) == "Rcpp_gbp3q")) {
    stop("gbp3q_viewer: sn must be a gbp3q object with class Rcpp_gbp3q.")
  }

  kk <- which(sn$k == 1)

  it <- data.table(
    id = kk,
    x  = sn$it[1, kk],
    y  = sn$it[2, kk],
    z  = sn$it[3, kk],
    l  = sn$it[4, kk],
    d  = sn$it[5, kk],
    h  = sn$it[6, kk]
  )

  ff <- which(sn$f == 1)

  bn <- data.table(
    id = ff,
    l  = sn$bn[1, ff],
    d  = sn$bn[2, ff],
    h  = sn$bn[3, ff]
  )

  if (is.null(title)) {
    title <- "gbp3q_viewer"
  }

  if (is.null(subtitle)) {
    subtitle <- paste0("bin: (", bn[["l"]], ", ", bn[["d"]], ", ", bn[["h"]], "): # fit: ", length(kk), " of ", length(sn$k))
  }

  bpp_viewer_single(
    it, bn, title = title, subtitle = subtitle
  )

}

#------------------------------------------------------------------------------#
