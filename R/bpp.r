#------------------------------------------------------------------------------#
#--------------------------------- gbp::gbp.r ---------------------------------#
#------------------------- author: gyang274@gmail.com -------------------------#
#------------------------------------------------------------------------------#

#--------+---------+---------+---------+---------+---------+---------+---------#
#234567890123456789012345678901234567890123456789012345678901234567890123456789#

#------------------------------------------------------------------------------#
#------------------------------------ main ------------------------------------#
#------------------------------------------------------------------------------#

#' gbp
#'
#' @description
#'
#'  a collection of 1d, 2d, 3d and 4d bin packing problem solver
#'
#' @section solver:
#'
#'  r-level:
#'
#'   wrapper over c-level function aims solving e-commerce bin packing problem
#'
#'    bpp_solver
#'
#'  c-level:
#'
#'   core class and solver on 1d, 2d, 3d and 4d bpp
#'
#'    gbp1d_solver_dpp
#'
#'    gbp2d_solver_dpp
#'
#'    gbp2d_solver_dpp_filt
#'
#'    gbp3d_solver_dpp
#'
#'    gbp3d_solver_dpp_filt
#'
#'    gbp4d_solver_dpp
#'
#'    gbp4d_solver_dpp_filt
#'
#'    bpp_solver_sgl
#'
#'    bpp_solver_dpp
#'
#' @section optimizer:
#'
#'  TODO: implementing a bin-shuffing optimizer?
#'
#' @section viewer:
#'
#'  rgl 3d show packing obtained via bpp_solver
#'
#'   bpp_viewer
#'
#' @docType package
#'
#' @name gbp
NULL

#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#------------------------------------ load ------------------------------------#
#------------------------------------------------------------------------------#
# loadModule: load on c++ codes class and module when package is loaded.
loadModule("gbp1d_module", TRUE)
loadModule("Ktlist2d_module", TRUE)
loadModule("gbp2d_module", TRUE)
loadModule("gbp2q_module", TRUE)
loadModule("Ktlist3d_module", TRUE)
loadModule("gbp3d_module", TRUE)
loadModule("gbp3q_module", TRUE)
loadModule("Ktlist4d_module", TRUE)
loadModule("gbp4d_module", TRUE)
loadModule("gbp4q_module", TRUE)
loadModule("bppSgl_module", TRUE)
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#----------------------------------- unload -----------------------------------#
#------------------------------------------------------------------------------#
# .onUnload: clean up c++ code used in package when package is unloaded.
.onUnload <- function (libpath) {
  library.dynam.unload("gbp", libpath)
}
#------------------------------------------------------------------------------#

