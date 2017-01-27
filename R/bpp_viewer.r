#------------------------------------------------------------------------------#
#--------------------------------- gbp::gbp.r ---------------------------------#
#------------------------- author: gyang274@gmail.com -------------------------#
#------------------------------------------------------------------------------#

#--------+---------+---------+---------+---------+---------+---------+---------#
#234567890123456789012345678901234567890123456789012345678901234567890123456789#

#------------------------------------------------------------------------------#
#------------------------------------ main ------------------------------------#
#------------------------------------------------------------------------------#

#' bpp_viewer
#'
#' @description
#'
#'  bpp single or multiple order packing solution viewer
#'
#' @param sn
#'
#'  sn bpp_solution from bpp_solver <list>
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
#' @param title
#'
#'  title <character>
#'
#' @param subtitle
#'
#'  subtitle <character>
#'
#' @family bpp_viewer
#' @rdname bpp_viewer
#' @export
bpp_viewer <- function(
  sn, title = NULL, subtitle = NULL
) {

  if (is.null(title)) {
    title <- "bpp_viewer"
  }

  it <- sn[["it"]]

  it <- it[ , list(
    otid = get("otid"),
    bid  = get("bid"),
    id   = get("sku"),
    x    = get("x"),
    y    = get("y"),
    z    = get("z"),
    l    = get("l"),
    d    = get("d"),
    h    = get("h"),
    w    = get("w")
  )]

  bn <- sn[["bn"]]

  lktd <- unique(it[["otid"]])

  for (i in lktd) {

    lkit <- it[get("otid") == i]

    lkbd <- unique(lkit[["bid"]])

    if (length(lkbd) != 1) {
      warning("bpp_viewer: should fit one ticket into one bin?\n"); break;
    }

    lkbn <- bn[id == lkbd]

    bpp_viewer_single(it = lkit, bn = lkbn, title = title, subtitle = subtitle)

  }

}

#' bpp_viewer_single
#'
#' @description
#'
#'  bpp solution viewer on single bin all item
#'
#' @param it
#'
#'  it item <data.table>
#'
#'  - id it id <integer>
#'
#'  - x, y, z it position w.r.t bins <numeric>
#'
#'  - l, d, h it scale along x, y, z <numeric>
#'
#'  - w it weight <numeric>
#'
#'  - auto: cc, wd, txt point and lines color, size, legend <numeric/character, numeric, character>
#'
#' @param bn
#'
#'  bn bins <data.table>
#'
#'  - id bn id <integer>
#'
#'  - l, d, h bn scale <numeric>
#'
#'  - w bn weight limit <numeric>
#'
#'  - auto: cc, wd, txt point and lines color, size, legend <numeric/character, numeric, character>
#'
#' @param title
#'
#'  title <character>
#'
#' @param subtitle
#'
#'  subtitle <character>
#'
#' @param it_rgl_control
#'
#'  control the color of it in rgl
#'
#' @param bn_rgl_control
#'
#'  control the color of bn in rgl
#'
#' @param label_it
#'
#'  label text on it corner or not
#'
#' @param label_bn
#'
#'  label text on bn corner or not
#'
#' @family bpp_viewer
#' @rdname bpp_viewer_single
#' @export
bpp_viewer_single <- function(
  it, bn, title = NULL, subtitle = NULL, it_rgl_control = NULL, bn_rgl_control = NULL, label_it = TRUE, label_bn = TRUE
) {

  #- init
  if (nrow(bn) != 1L) {
    stop("bpp_viewer_single: bn must have 1 row.\n")
  }

  if (is.null(it_rgl_control)) {
    it_rgl_control <- create_it_rgl_control()
  }

  if (is.null(bn_rgl_control)) {
    bn_rgl_control <- create_bn_rgl_control()
  }

  #- auto: add color size and legend to it
  nm <- names(it)

  if (!("id" %in% nm)) {
    it[ , `:=`(id = 1:.N)]
  }

  if (!("cc" %in% nm)) {
    it[ , `:=`(cc = it_rgl_control[["cc"]][c(c(1L:.N - 1L) %% length(it_rgl_control[["cc"]]) + 1L)])]
  }

  if (!("wd" %in% nm)) {
    it[ , `:=`(wd = it_rgl_control[["wd"]][c(c(1L:.N - 1L) %% length(it_rgl_control[["wd"]]) + 1L)])]
  }

  if (!("txt" %in% nm)) {
    if (!("w" %in% nm)) {
      it[ , `:=`(txt = paste0("item-", id))]
    } else {
      it[ , `:=`(txt = paste0("item-", id, " (W: ", round(w, digits = 2L), ")"))]
    }
  }

  #- auto: add color size and legend to bn
  mn <- names(bn)

  if (!("id" %in% mn)) {
    bn[ , `:=`(id = 1:.N)]
  }

  if (!(all(c("x", "y", "z") %in% mn))) {
    bn[ , `:=`(x = 0, y = 0, z = 0)]
  }

  if (!("cc" %in% mn)) {
    bn[ , `:=`(cc = bn_rgl_control[["cc"]][c(c(1L:.N - 1L) %% length(bn_rgl_control[["cc"]]) + 1L)])]
  }

  if (!("wd" %in% mn)) {
    bn[ , `:=`(wd = bn_rgl_control[["wd"]][c(c(1L:.N - 1L) %% length(bn_rgl_control[["wd"]]) + 1L)])]
  }

  if (!("txt" %in% mn)) {
    bn[ , `:=`(txt = paste0("bin-", id))]
  }

  #- main
  rgl::open3d()

  # rgl::bg3d("#f0f2f4")

  rgl::par3d("windowRect" = c(120, 120, 600, 600))

  #- it cube3d
  for (i in 1L:nrow(it)) {
    create_it_cube3d(
      it[i, get("id")], it[i, get("x")], it[i, get("y")], it[i, get("z")], it[i, get("l")], it[i, get("d")], it[i, get("h")],
      it[i, get("cc")], it[i, get("wd")], it[i, get("txt")], TRUE
    )
  }

  #- bn cube3d
  create_it_cube3d(
    bn[1L, get("id")], bn[1L, get("x")], bn[1L, get("y")], bn[1L, get("z")], bn[1L, get("l")], bn[1L, get("d")], bn[1L, get("h")],
    bn[1L, get("cc")], bn[1L, get("wd")], bn[1L, get("txt")], FALSE
  )

  #- add background
  rgl::bgplot3d({

    graphics::par(bg = "#f0f2f4")

    graphics::plot.new()
    graphics::title(main = ifelse(is.null(title), "bpp_viewer_single", title), line = 3)
    graphics::mtext(side = 1, ifelse(is.null(subtitle), paste0("bin: ", bn[["id"]], " (", bn[["l"]], ", ", bn[["d"]], ", ", bn[["h"]], "): ", ifelse("w" %in% nm, paste0(", total weight: ", sum(it[["w"]])), "")), subtitle), line = 4)

    # graphics::par(bg = "#ffffff") # set back bg?
  })

}

#' create_it_cube3d
#' @description
#'  subroutine of bpp_viewer_single
#' @details
#'  add it or bn on current rgl device
#' @param id id
#' @param x x-coordinate
#' @param y y-coordinate
#' @param z z-coordinate
#' @param l length along x-coordinate
#' @param d depth  along y-coordinate
#' @param h height along z-coordinate
#' @param cc color
#' @param wd width
#' @param txt text
#' @param itxt plot text or not
#' @export
create_it_cube3d <- function(id, x, y, z, l, d, h, cc, wd, txt, itxt = TRUE) {

  #- it cube3d vertex via coordinate
  it <- data.table::data.table(
    id = rep.int(id, times = 8L),
    x  = rep.int(x, times = 8L) + rep(c(0, l), times = 4L),
    y  = rep.int(y, times = 8L) + rep(rep(c(0, d), each = 2L), times = 2L),
    z  = rep.int(z, times = 8L) + rep(c(0, h), each = 4L),
    cc = rep.int(cc, times = 8L)
  )

  # vectorize with single size
  it %$% rgl::points3d(x, y, z, col = cc, size = wd)

  #- it cube3d edge via line segment

  #- mapping 1 point point3d in id -> segment3d in jd

  id <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L)

  jd <- c(1L, 2L, 1L, 3L, 1L, 5L, 2L, 4L, 2L, 6L, 3L, 4L, 3L, 7L, 4L, 8L, 5L, 6L, 5L, 7L, 6L, 8L, 7L, 8L)

  it[jd, ] %$% rgl::segments3d(x, y, z, col = cc, size = wd)

  #- it cube3d label via text on xyz
  if (itxt) {
    rgl::text3d(x - 0.50, y - 0.50, z - 0.50, texts = txt, col = cc)
  }

}

#' create_it_rgl_control
#' @description
#'  subroutine of bpp_viewer_single
#' @export
create_it_rgl_control <- function() {

  list(
    cc = c("#dd4b39", "#ffaf16", "#f39c12", "#00a65a", "#00c0ef", "#337ab7", "#605ca8"), wd = 8L
  )

}

#' create_bn_rgl_control
#' @description
#'  subroutine of bpp_viewer_single
#' @export
create_bn_rgl_control <- function() {

  list(
    # cc = c("#d2d6de"), wd = 8L
    # cc = c("#7b879d"), wd = 8L
    cc = c("#21242c"), wd = 8L
  )

}
#------------------------------------------------------------------------------#

