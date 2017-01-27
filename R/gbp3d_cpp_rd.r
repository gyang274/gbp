#' Ktlist3d
#' @aliases
#'  Ktlist3d Rcpp_Ktlist3d Rcpp_Ktlist3d-class
#' @description
#'  Ktlist3d hold multiple kt for recursive fit
#' @details
#'  Ktlist3d hold multiple kt via consider all possible fit onto different xp and different rotation and nlimit
#'
#'  a Ktlist3d class instance has 4 fields:
#'
#'   - n: length of kt candidate position scale vector list
#'
#'   - kt: candidate (x, y, z, l, d, h) fit of it current investigating
#'
#'   - xp: candidate extreme point list after kt fit into each corresponding (x, y, z, l, d, h) position scale
#'
#'   - s: score of each kt fit: calculate overall extrem point residual space entropy as score, the smaller the better,
#'   since smaller entropy indicate concentrated residual space and less number of extreme point.
#' @note
#'  internal cpp class use in gbp3d_solver_dpp()
#' @family gbp3d_it
#' @rdname Ktlist3d
#' @docType class
"Ktlist3d"


#' gbp3d_it_create_ktlist
#' @description
#'  create ktlist from itlist
#' @details
#'  core function in gbp3d_solver_dpp
#'   select highest profitable it not yet fit into bn and return all possbile fit w.r.t xp and orientation
#' @param bn
#'  bn scale <vector>
#'  - l, d, h bn scale along x, y, z <numeric>
#' @param it
#'  it position and scale <matrix>
#'  - x, y, z it position in the bin <numeric>
#'  - l, d, h it scale along x, y, z <numeric>
#' @param xp
#'  xp extreme point position and residual space scale <matrix>
#'  - x, y, z xp position in the bin <numeric>
#'  - l, d, h xp residual space scale along x, y, z <numeric>
#' @param ktinit
#'  kt candidate scale without position <matrix>
#'  - l, d, h kt scale along x, y, z which open to orientation <numeric>
#' @param nlmt
#'  nlmt: limit on ktlist n max-value
#' @return Ktlist3d
#' @note
#'  should make sure it kt can be fit in bin outside
#' @note
#'  internal function use in gbp3d_solver_dpp() for creating Ktlist3d object for fit.
#' @family gbp3d_it
#' @rdname gbp3d_it_create_ktlist
"gbp3d_it_create_ktlist"


#' gbp3d
#' @aliases
#'  gbp3d Rcpp_gbp3d Rcpp_gbp3d-class
#' @description
#'  generalized bin packing problem in 3 dimension, a.k.a bin packing problem.
#' @details
#'  gbp3d init a profit vector p, a length vector l, a depth vector d, a height vector h, and also
#'   a length constraint ml, a depth constraint md, and a height constraint mh on l x d x h cuboid
#'   with geometry intepretation.
#'
#'  gbp3d solver would solve
#'
#'    maximize   sum_{j=1}^{n} p_{j} k_{j}
#'
#'    subject to fit (l_{j}, d_{j}, h_{j}) at coordinate (x_{j}, y_{j}, z_{j})
#'               such that no overlap in ml x md x mh cuboid, j = 1, ......, n
#'
#'  and instantiate a gbp3d object with a x-axis coordinate vector x, a y-axis coordinate vector y,
#'   a z-axis coordinate vector z, a selection vector k, and an objective o.
#'
#'  a gbp3d class instance has 6 fields:
#'
#'   - p: profit of it fit into bn <vector>
#'
#'        created via cluster max(l, d, h) and area via gbp3d_solver_dpp_main_create_p()
#'
#'   - it: it position and scale <matrix>
#'
#'     - x, y, z it position in the bin <numeric>
#'
#'     - l, d, h it scale along x, y, z <numeric>
#'
#'   - bn: bn scale <vector>
#'
#'     - l, d, h bn scale along x, y, z <numeric>
#'
#'   - k: selection indicator 0, 1 <vector>
#'
#'   - o: objective achivement volumn fit in over volumn overall <numeric>
#'
#'   - ok: a quick indicator of all it fit into bn? <bool>
#'
#' @note
#'  p is a proxy of ranking on cuboid fit difficulty, often a func of max(l, d, h), surface, volume
#'   and solver would often maximize sum_{j=1}^{n} v_{j} k_{j} instead of sum_{j=1}^{n} p_{j} k_{j}
#' @family gbp3d
#' @rdname gbp3d
#' @docType class
#' @export
"gbp3d"


#' gbp3d_solver_dpp
#' @description
#'  solve gbp3d via extreme point heuristic and best information score fit strategy.
#' @details
#'  gbp3d init a profit vector p, a length vector l, a depth vector d, a height vector h, and also
#'   a length constraint ml, a depth constraint md, and a height constraint mh on l x d x h cuboid
#'   with geometry intepretation.
#'
#'  gbp3d solver would solve
#'
#'    maximize   sum_{j=1}^{n} p_{j} k_{j}
#'
#'    subject to fit (l_{j}, d_{j}, h_{j}) at coordinate (x_{j}, y_{j}, z_{j})
#'               such that no overlap in ml x md x mh cuboid, j = 1, ......, n
#'
#'  and instantiate a gbp3d object with a x-axis coordinate vector x, a y-axis coordinate vector y,
#'   a z-axis coordinate vector z, a selection vector k, and an objective o.
#' @param p
#'  p profit of it fit into bn <vector>
#'  - cluster max(l, d) and min(l, d) via gbp3d_solver_dpp_prep_create_p()
#' @param ldh
#'  it position and scale <matrix>
#'  - l, d, h it scale along x, y, z, subject to orientation rotation <numeric>
#' @param m
#'  bn scale <vector>
#'  - l, d, h bn scale along x, y, z <numeric>
#' @return gbp3d
#'  a gbp3d instantiate with p profit, it item (x, y, z, l, d, h) position scale matrix, bn bin (l, d, h) scale vector,
#'   k selection, o objective, and ok an indicator of all fit or not.
#' @family gbp3d
#' @rdname gbp3d_solver_dpp
"gbp3d_solver_dpp"


#' gbp3d_checkr
#' @description
#'  auxilium of gbp3d_solver_dpp
#' @details
#'  check fit solution is valid:
#'   no conflict between item and bin, and no conflict between each pair of item.
#' @param sn <gbp3d>
#'  gbp3d object from gbp3d_solver_dpp() solution.
#' @return okfit? <bool>
#' @family gbp3d
#' @rdname gbp3d_checkr
"gbp3d_checkr"


#' gbp3q
#' @aliases
#'  gbp3q Rcpp_gbp3q Rcpp_gbp3q-class
#' @description
#'  generalized bin packing problem in 3 dimension, a.k.a bin packing problem.
#' @details
#'  gbp3d init a profit vector p, a length vector l, a depth vector d, a height vector h, and also
#'   a length constraint ml, a depth constraint md, and a height constraint mh on l x d x h cuboid
#'   with geometry intepretation.
#'
#'  gbp3d solver would solve
#'
#'    maximize   sum_{j=1}^{n} p_{j} k_{j}
#'
#'    subject to fit (l_{j}, d_{j}, h_{j}) at coordinate (x_{j}, y_{j}, z_{j})
#'               such that no overlap in ml x md x mh cuboid, j = 1, ......, n
#'
#'  and instantiate a gbp3d object with a x-axis coordinate vector x, a y-axis coordinate vector y,
#'   a z-axis coordinate vector z, a selection vector k, and an objective o.
#'
#'  gbp3q solver would also select the most preferred often smallest m from a list of m(l, d, h) after
#'   determine all or the higest volume set of ld can fit into one m(l, d, h).
#'
#'  a gbp3q class instance has 7 fields:
#'
#'   - p: profit of it fit into bn <vector>
#'
#'        created via cluster max(l, d, h) and area via gbp3d_solver_dpp_main_create_p()
#'
#'   - it: it position and scale <matrix>
#'
#'     - x, y, z it position in the bin <numeric>
#'
#'     - l, d, h it scale along x, y, z <numeric>
#'
#'   - bn: bn scale <matrix>
#'
#'     - l, d, h bn scale along x, y, z <numeric>
#'
#'       matrix of 3 rows and each column is a single bn
#'
#'     should make sure bn list are sorted via volume
#'      so that the first col is the most prefered smallest bn, and also
#'      the last col is the least prefered largest and often dominant bn
#'
#'     should make sure no X in front of Y if bnX dominant bnY,
#'      bnX dominant bnY if all(X(l, d, h) > Y(l, d, h)) and should always prefer Y.
#'
#'     should make sure bn such that l >= d or vice versa.
#'
#'   - k: selection indicator 0, 1 on it <vector>
#'
#'   - f: selection indicator 0, 1, 2, 3 on bn <vector>
#'
#'        f in result should have no 0 and only one of 1
#'
#'   - o: objective achivement volumn fit in over volumn overall <numeric>
#'
#'   - ok: a quick indicator of all it fit into bn? <bool>
#'
#' @family gbp3q
#' @rdname gbp3q
#' @docType class
#' @export
"gbp3q"


#' gbp3d_solver_dpp_filt
#' @description
#'  solve gbp3d w.r.t select most preferable often smallest bin from bn list
#' @details
#'  gbp3d_solver_dpp_filt is built on top of gbp3d_solver_dpp
#'   aims to select the most preferable bn from a list of bn that can fit all or most it
#'
#'  gbp3d_solver_dpp()'s objective is fit all or most it into a single given bn (l, d, h)
#'
#'  gbp3d_solver_dpp_filt()'s objective is select the most preferable given a list of bn
#'   where bn list is specified in 3xN matrix that the earlier column the more preferable
#'
#'  gbp3d_solver_dpp_filt() use an approx binary search and determine f w.r.t bn.n_cols
#'   where f = 1 indicate the bn being selected and only one of 1 in result returned.
#'
#'  ok = true if any bin can fit all it and algorithm will select smallest bn can fit all
#'   otherwise ok = false and algorithm will select a bn can maximize volume of fitted it
#'
#'  often recommend to make the last and least preferable bn dominate all other bn in list
#'   when design bn list, bnX dominant bnY if all(X(l, d, h) > Y(l, d, h)).
#' @param ldh
#'  it scale <matrix>
#'  - l, d, h it scale along x, y, z <numeric>
#' @param m
#'  bn scale <matrix>
#'  - l, d, h bn scale along x, y, z <numeric>
#'  - l, d, h in row and each col is a single bn
#'
#'  should make sure bn list are sorted via volume
#'   so that the first col is the most prefered smallest bn, and also
#'   the last col is the least prefered largest and often dominant bn
#'
#'  should make sure no X in front of Y if bnX dominant bnY,
#'   bnX dominant bnY if all(X(l, d, h) > Y(l, d, h)) and should always prefer Y.
#'
#'  should make sure bn such that l >= d >= h or vice versa.
#' @return gbp3q
#'  a gbp3q instantiate with p profit, it item (x, y, z, l, d, h) position scale matrix, bn bin (l, d, h) scale matrix,
#'   k it selection, o objective, f bn selection, and ok an indicator of all fit or not.
#' @family gbp3q
#' @rdname gbp3d_solver_dpp_filt
"gbp3d_solver_dpp_filt"


#' gbp3q_checkr
#' @description
#'  auxilium of gbp3d_solver_dpp_filt
#' @details
#'  check fit solution is valid:
#'   no conflict between item and bin, and no conflict between each pair of item.
#' @param sn <gbp3q>
#'  gbp3q object from gbp3d_solver_dpp_filt() solution.
#' @return okfit? <bool>
#' @family gbp3q
#' @rdname gbp3q_checkr
"gbp3q_checkr"

