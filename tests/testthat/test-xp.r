# testthat::context("extreme point calculation")

# testthat::test_that(
#   "xp: extreme point calculation - projecting new extreme point to it surface", {
#     testthat::expect_equal(
#       gbp3d_xp_it_pjt_kt(
#         it = c(0, 0, 0, 4, 4, 4),
#         kt = c(7, 0, 0, 2, 2, 2)
#       ),
#       matrix(c(0, 0, 0, 1, 1, 0), 6, 1)
#     )
#     testthat::expect_equal(
#       gbp3d_xp_it_pjt_kt(
#         it = c(0, 0, 0, 4, 4, 4),
#         kt = c(0, 7, 0, 2, 2, 2)
#       ),
#       matrix(c(1, 0, 0, 0, 0, 1), 6, 1)
#     )
#     testthat::expect_equal(
#       gbp3d_xp_it_pjt_kt(
#         it = c(0, 0, 0, 4, 4, 4),
#         kt = c(0, 0, 7, 2, 2, 2)
#       ),
#       matrix(c(0, 1, 1, 0, 0, 0), 6, 1)
#     )
#   }
# )
