#' Ktlist4d
#' @aliases
#'  Ktlist4d Rcpp_Ktlist4d Rcpp_Ktlist4d-class
#' @description
#'  Ktlist4d hold multiple kt for recursive fit
#' @details
#'  Ktlist4d hold multiple kt via consider all possible fit onto different xp and different rotation and nlimit
#'
#'  a Ktlist4d class instance has 4 fields
#'
#'   - n: length of kt candidate position scale vector list
#'
#'   - kt: candidate (x, y, z, w, l, d, h, w) fit of it current investigating
#'
#'         x, y, z, w - weight holding in bn when fit; l, d, h, w - weight of it itself
#'
#'   - xp: candidate extreme point list after kt fit into each corresponding (x, y, z, l, d, h) position scale
#'
#'   - s: score of each kt fit: calculate overall extrem point residual space entropy as score, the smaller the better,
#'   since smaller entropy indicate concentrated residual space and less number of extreme point.
#' @note
#'  internal cpp class use in gbp4d_solver_dpp()
#' @family gbp4d_it
#' @rdname Ktlist4d
#' @docType class
"Ktlist4d"


#' gbp4d_it_create_ktlist
#' @description
#'  create ktlist from itlist
#' @details
#'  core function in gbp4d_solver_dpp
#'   select highest profitable it not yet fit into bn and return all possbile fit w.r.t xp and orientation
#' @param bn
#'  bn scale <vector>
#'  - l, d, h, w bn scale along x, y, z and w <numeric>
#' @param it
#'  it position and scale <matrix>
#'  - x, y, z, w it position and w in the bin <numeric>
#'  - l, d, h, w it scale along x, y, z and w <numeric>
#' @param xp
#'  xp extreme point position and residual space scale <matrix>
#'  - x, y, z, w xp position and w in the bin <numeric>
#'  - l, d, h, w xp residual space scale along x, y, z and w <numeric>
#' @param ktinit
#'  kt candidate scale without position <matrix>
#'  - l, d, h, w kt scale along x, y, z, w which open to orientation <numeric>
#' @param nlmt
#'  nlmt: limit on ktlist n max-value
#' @return Ktlist4d
#' @note
#'  should make sure it kt can be fit in bin outside
#' @note
#'  internal function use in gbp4d_solver_dpp() for creating Ktlist4d object for fit.
#' @family gbp4d_it
#' @rdname gbp4d_it_create_ktlist
"gbp4d_it_create_ktlist"


#' gbp4d
#' @aliases
#'  gbp4d Rcpp_gbp4d Rcpp_gbp4d-class
#' @description
#'  generalized bin packing problem in 4 dimension, a.k.a bin packing problem with weight limit.
#' @details
#'  gbp4d init a profit vector p, a length l, a depth d, a height h, and a weight w, along with
#'   associate constraints ml, md, mh and mw.
#'  gbp4d should fit it (l, d, h, w) into bn (ml, md, mh, mw) with w on weight limit constraint
#'   and l, d, h on geometry intepretation.
#'  gbp4d solver would solve
#'
#'    maximize   sum_{j=1}^{n} p_{j} k_{j}
#'
#'    subject to sum_{j=1}^{n} w_{j} k_{j} leq mw and
#'
#'               fit (l_{j}, d_{j}, h_{j}) at coordinate (x_{j}, y_{j}, z_{j})
#'               such that no overlap in ml x md x mh cuboid, j = 1, ......, n
#'
#'  and instantiate a gbp4d object with a x-axis coordinate vector x, a y-axis coordinate vector y,
#'   a z-axis coordinate vector z, a selection vector k, and an objective o.
#'
#'  a gbp4d class instance has 6 fields:
#'
#'   - p: profit of it fit into bn <vector>
#'
#'        created via cluster w via gbp1d, cluster max(l, d, h) and area via gbp4d_solver_dpp_main_create_p()
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
#'   - k: selection indicator 0, 1 <vector>
#'
#'   - o: objective achivement volumn fit in over volumn overall <numeric>
#'
#'   - ok: a quick indicator of all it fit into bn? <bool>
#'
#' @note
#'  p is a proxy of ranking on cuboid fit difficulty, often a func of max(l, d, h), surface, volume
#'   and solver would often maximize sum_{j=1}^{n} v_{j} k_{j} instead of sum_{j=1}^{n} p_{j} k_{j}
#' @family gbp4d
#' @rdname gbp4d
#' @docType class
#' @export
"gbp4d"


#' gbp4d_solver_dpp
#' @description
#'  solve gbp4d via extreme point heuristic and best information score fit strategy.
#' @details
#'  gbp4d init a profit vector p, a length l, a depth d, a height h, and a weight w, along with
#'   associate constraints ml, md, mh and mw.
#'  gbp4d should fit it (l, d, h, w) into bn (ml, md, mh, mw) with w on weight limit constraint
#'   and l, d, h on geometry intepretation.
#'  gbp4d solver would solve
#'
#'    maximize   sum_{j=1}^{n} p_{j} k_{j}
#'
#'    subject to sum_{j=1}^{n} w_{j} k_{j} leq mw and
#'
#'               fit (l_{j}, d_{j}, h_{j}) at coordinate (x_{j}, y_{j}, z_{j})
#'               such that no overlap in ml x md x mh cuboid, j = 1, ......, n
#'
#'  and instantiate a gbp4d object with a x-axis coordinate vector x, a y-axis coordinate vector y,
#'   a z-axis coordinate vector z, a selection vector k, and an objective o.
#' @param p
#'  p profit of it fit into bn <vector>
#'  - cluster w via gbp1d, cluster max(l,d,h) and area via gbp4d_solver_dpp_main_create_p()
#' @param ldhw
#'  it scales <matrix>
#'  - l, d, h, w it scale along x, y, z and w (weight on separate single dimension) <numeric>
#' @param m
#'  bn scales <vector>
#'  - l, d, h, w bn scale along x, y, z and w (weight on separate single dimension) <numeric>
#' @return gbp4d
#'  a gbp4d instantiate with p profit, it item (x, y, z, w, l, d, h, w) position scale matrix, bn bin (l, d, h, w) scale vector,
#'   k selection, o objective, and ok an indicator of all fit or not.
#' @family gbp4d
#' @rdname gbp4d_solver_dpp
"gbp4d_solver_dpp"


#' gbp4d_checkr
#' @description
#'  auxilium of gbp4d_solver_dpp
#' @details
#'  check fit solution is valid:
#'   no conflict between item and bin, and no conflict between each pair of item, and no conflict on weight limit.
#' @param sn <gbp4d>
#'  gbp4d object from gbp4d_solver_dpp() solution.
#' @return okfit? <bool>
#' @family gbp4d
#' @rdname gbp4d_checkr
"gbp4d_checkr"


#' gbp4q
#' @aliases
#'  gbp4q Rcpp_gbp4q Rcpp_gbp4q-class
#' @description
#'  generalized bin packing problem in 4 dimension, a.k.a bin packing problem with weight limit.
#' @details
#'  gbp4d init a profit vector p, a length l, a depth d, a height h, and a weight w, along with
#'   associate constraints ml, md, mh and mw.
#'  gbp4d should fit it (l, d, h, w) into bn (ml, md, mh, mw) with w on weight limit constraint
#'   and l, d, h on geometry intepretation.
#'  gbp4d solver would solve
#'
#'    maximize   sum_{j=1}^{n} p_{j} k_{j}
#'
#'    subject to sum_{j=1}^{n} w_{j} k_{j} leq mw and
#'
#'               fit (l_{j}, d_{j}, h_{j}) at coordinate (x_{j}, y_{j}, z_{j})
#'               such that no overlap in ml x md x mh cuboid, j = 1, ......, n
#'
#'  and instantiate a gbp4d object with a x-axis coordinate vector x, a y-axis coordinate vector y,
#'   a z-axis coordinate vector z, a selection vector k, and an objective o.
#'
#'  gbp4q solver would also select the most preferred often smallest m from a list of m(l, d, h) after
#'   determine all or the higest volume set of ld can fit into one m(l, d, h) w.r.t the weight constraint.
#'
#'  a gbp4q class instance has 7 fields:
#'
#'   - p: profit of it fit into bn <vector>
#'
#'        created via cluster w via gbp1d, cluster max(l, d, h) and area via gbp4d_solver_dpp_main_create_p()
#'
#'   - it: it position and scale <matrix>
#'
#'     - x, y, z, w it position and w in the bin <numeric> (w hold in bn when fit it in bn)
#'
#'     - l, d, h, w it scale along x, y, z and w <numeric> (w of it itself)
#'
#'   - bn: bn scale <matrix>
#'
#'     - l, d, h, w bn scale along x, y, z and w <numeric>
#'
#'       matrix of 4 rows and each column is a single bn
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
#' @family gbp4q
#' @rdname gbp4q
#' @docType class
#' @export
"gbp4q"


#' gbp4d_solver_dpp_filt
#' @description
#'  solve gbp4d w.r.t select most preferable often smallest bin from bn list
#' @details
#'  gbp4d_solver_dpp_filt is built on top of gbp4d_solver_dpp
#'   aims to select the most preferable bn from a list of bn that can fit all or most it
#'
#'  gbp4d_solver_dpp()'s objective is fit all or most it into a single given bn (l, d, h, w)
#'
#'  gbp4d_solver_dpp_filt()'s objective is select the most preferable given a list of bn
#'   where bn list is specified in 4xN matrix that the earlier column the more preferable
#'
#'  gbp4d_solver_dpp_filt() use an approx binary search and determine f w.r.t bn.n_cols
#'   where f = 1 indicate the bn being selected and only one of 1 in result returned.
#'
#'  ok = true if any bin can fit all it and algorithm will select smallest bn can fit all
#'   otherwise ok = false and algorithm will select a bn can maximize volume of fitted it
#'
#'  often recommend to make the last and least preferable bn dominate all other bn in list
#'   when design bn list, bnX dominant bnY if all(X(l, d, h) > Y(l, d, h)).
#' @param ldhw
#'  it scale <matrix>
#'  - l, d, h, w it scale along x, y, z and w <numeric>
#' @param m
#'  bn scale <matrix>
#'  - l, d, h, w bn scale along x, y, z and w <numeric>
#'  - l, d, h, w in row and each col is a single bn
#'
#'  should make sure bn list are sorted via volume
#'   so that the first col is the most prefered smallest bn, and also
#'   the last col is the least prefered largest and often dominant bn
#'
#'  should make sure no X in front of Y if bnX dominant bnY,
#'   bnX dominant bnY if all(X(l, d, h) > Y(l, d, h)) and should always prefer Y.
#'
#'  should make sure bn such that l >= d >= h or vice versa.
#' @return gbp4q
#'  a gbp4q instantiate with p profit, it item (x, y, z, w, l, d, h, w) position scale matrix, bn bin (l, d, h, w) scale matrix,
#'   k it selection, o objective, f bn selection, and ok an indicator of all fit or not.
#' @family gbp4q
#' @rdname gbp4d_solver_dpp_filt
"gbp4d_solver_dpp_filt"


#' gbp4q_checkr
#' @description
#'  auxilium of gbp4d_solver_dpp_filt
#' @details
#'  check fit solution is valid:
#'   no conflict between item and bin, and no conflict between each pair of item, and no conflict on weight limit.
#' @param sn <gbp4q>
#'  gbp4q object from gbp4d_solver_dpp_filt() solution.
#' @return okfit? <bool>
#' @family gbp4q
#' @rdname gbp4q_checkr
"gbp4q_checkr"

