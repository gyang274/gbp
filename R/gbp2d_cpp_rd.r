#' Ktlist2d
#' @aliases
#'  Ktlist2d Rcpp_Ktlist2d Rcpp_Ktlist2d-class
#' @description
#'  Ktlist2d hold multiple kt for recursive fit
#' @details
#'  Ktlist2d hold multiple kt via consider all possible fit onto different xp and different rotation and nlimit
#'
#'  a Ktlist2d class instance has 4 fields:
#'
#'   - n: length of kt candidate position scale vector list
#'
#'   - kt: candidate (x, y, l, d) fit of it current investigating
#'
#'   - xp: candidate extreme point list after kt fit into each corresponding (x, y, l, d) position scale
#'
#'   - s: score of each kt fit: calculate overall extrem point residual space entropy as score, the smaller the better,
#'   since smaller entropy indicate concentrated residual space and less number of extreme point.
#' @note
#'  internal cpp class use in gbp2d_solver_dpp()
#' @family gbp2d_it
#' @rdname Ktlist2d
#' @docType class
#' @export
"Ktlist2d"


#' gbp2d_it_create_ktlist
#' @description
#'  create ktlist from itlist
#' @details
#'  core function in gbp2d_solver_dpp
#'   select highest profitable it not yet fit into bn and return all possbile fit w.r.t xp and orientation
#' @param bn
#'  bn scale <vector>
#'  - l, d bn scale along x and y <numeric>
#' @param it
#'  it position and scale <matrix>
#'  - x, y it position in the bin <numeric>
#'  - l, d it scale along x and y <numeric>
#' @param xp
#'  xp extreme point position and residual space scale <matrix>
#'  - x, y xp position in the bin <numeric>
#'  - l, d xp residual space scale along x and y <numeric>
#' @param ktinit
#'  kt candidate scale without position <matrix>
#'  - l, d kt scale along x and y which open to orientation <numeric>
#' @param nlmt
#'  nlmt: limit on ktlist n max-value
#' @return Ktlist2d
#' @note
#'  should make sure it kt can be fit in bin outside
#' @note
#'  internal function use in gbp2d_solver_dpp() for creating Ktlist2d object for fit.
#' @family gbp2d_it
#' @rdname gbp2d_it_create_ktlist
"gbp2d_it_create_ktlist"


#' gbp2d
#' @aliases
#'  gbp2d Rcpp_gbp2d Rcpp_gbp2d-class
#' @description
#'  generalized bin packing problem in 2 dimension, a.k.a rectangle fill.
#' @details
#'  gbp2d init a profit vector p, a length vector l, a depth vector d,
#'   a length constraint ml, and a depth constraint md on l x d rectangle
#'   with geometry intepretation.
#'
#'  gbp2d solver would solve
#'
#'    maximize   sum_{j=1}^{n} p_{j} k_{j}
#'
#'    subject to fit (l_{j}, d_{j}) at coordinate (x_{j}, y_{j})
#'               such that no overlap in ml x md, j = 1, ...., n
#'
#'  and instantiate a gbp2d object with a x-axis coordinate vector x, a y-axis coordinate vector y,
#'   a selection vector k, and an objective o.
#'
#'  a gbp2d class instance has 6 fields:
#'
#'   - p: profit of it fit into bn <vector>
#'
#'        created via cluster max(l, d) and min(l, d) via gbp2d_solver_dpp_prep_create_p()
#'
#'   - it: it position and scale <matrix>
#'
#'     - x, y it position in the bin <numeric>
#'
#'     - l, d it scale along x and y <numeric>
#'
#'   - bn: bn scale <vector>
#'
#'     - l, d bn scale along x and y <numeric>
#'
#'   - k: selection indicator 0, 1 <vector>
#'
#'   - o: objective achivement volumn fit in over volumn overall <numeric>
#'
#'   - ok: a quick indicator of all it fit into bn? <bool>
#'
#' @note
#'  p is a proxy of ranking on rectangle fit difficulty, often a function w.r.t max(l, d) and l x d
#' @family gbp2d
#' @rdname gbp2d
#' @docType class
#' @export
"gbp2d"


#' gbp2d_solver_dpp
#' @description
#'  solve gbp2d via extreme point heuristic and best information score fit strategy.
#' @details
#'  gbp2d init a profit vector p, a length vector l, a depth vector d,
#'   a length constraint ml, and a depth constraint md on l x d rectangle
#'   with geometry intepretation.
#'
#'  gbp2d solver would solve
#'
#'    maximize   sum_{j=1}^{n} p_{j} k_{j}
#'
#'    subject to fit (l_{j}, d_{j}) at coordinate (x_{j}, y_{j})
#'               such that no overlap in ml x md, j = 1, ...., n
#'
#'  and instantiate a gbp2d object with a x-axis coordinate vector x, a y-axis coordinate vector y,
#'   a selection vector k, and an objective o.
#' @param p
#'  p profit of it fit into bn <vector>
#'  - cluster max(l, d) and min(l, d) via gbp2d_solver_dpp_prep_create_p()
#' @param ld
#'  it position and scale <matrix>
#'  - l, d it scale along x and y, subject to orientation rotation <numeric>
#' @param m
#'  bn scale <vector>
#'  - l, d bn scale along x and y <numeric>
#' @return gbp2d
#'  a gbp2d instantiate with p profit, it item (x, y, l, d) position scale matrix, bn bin (l, d) scale vector,
#'   k selection, o objective, and ok an indicator of all fit or not.
#' @family gbp2d
#' @rdname gbp2d_solver_dpp
"gbp2d_solver_dpp"


#' gbp2d_checkr
#' @description
#'  auxilium of gbp2d and gbp2d_solver_dpp
#' @details
#'  check fit solution is valid:
#'   no conflict between item and bin, and no conflict between each pair of item.
#' @param sn <gbp2d>
#'  gbp2d object from gbp2d_solver_dpp() solution.
#' @return okfit? <bool>
#' @family gbp2d
#' @rdname gbp2d_checkr
"gbp2d_checkr"


#' gbp2q
#' @aliases
#'  gbp2q Rcpp_gbp2q Rcpp_gbp2q-class
#' @description
#'  generalized bin packing problem in 2 dimension, a.k.a rectangle fill.
#' @details
#'  gbp2d init a profit vector p, a length vector l, a depth vector d,
#'   a length constraint ml, and a depth constraint md on l x d rectangle
#'   with geometry intepretation.
#'
#'  gbp2d solver would solve
#'
#'    maximize   sum_{j=1}^{n} p_{j} k_{j}
#'
#'    subject to fit (l_{j}, d_{j}) at coordinate (x_{j}, y_{j})
#'               such that no overlap in ml x md, j = 1, ...., n
#'
#'  and instantiate a gbp2d object with a x-axis coordinate vector x, a y-axis coordinate vector y,
#'   a selection vector k, and an objective o.
#'
#'  gbp2q solver would also select the most preferred often smallest m from a list of m(l, d) after
#'   determine all or the higest volume set of ld can fit into one m(l, d).
#'
#'  a gbp2q class instance has 7 fields:
#'
#'   - p: profit of it fit into bn <vector>
#'
#'        created via cluster max(l, d) and min(l, d) via gbp2d_solver_dpp_prep_create_p()
#'
#'   - it: it position and scale <matrix>
#'
#'     - x, y it position in the bin <numeric>
#'
#'     - l, d it scale along x and y <numeric>
#'
#'   - bn: bn scale <matrix>
#'
#'     - l, d bn scale along x and y <numeric>
#'
#'       matrix of 2 rows and each column is a single bn
#'
#'     should make sure bn list are sorted via volume
#'      so that the first col is the most prefered smallest bn, and also
#'      the last col is the least prefered largest and often dominant bn
#'
#'     should make sure no X in front of Y if bnX dominant bnY,
#'      bnX dominant bnY if all(X(l, d) > Y(l, d)) and should always prefer Y.
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
#' @family gbp2q
#' @rdname gbp2q
#' @docType class
#' @export
"gbp2q"


#' gbp2d_solver_dpp_filt
#' @description
#'  solve gbp2d w.r.t select most preferable often smallest bin from bn list
#' @details
#'  gbp2d_solver_dpp_filt is built on top of gbp2d_solver_dpp
#'   aims to select the most preferable bn from a list of bn that can fit all or most it
#'
#'  gbp2d_solver_dpp()'s objective is fit all or most it into a single given bn (l, d)
#'
#'  gbp2d_solver_dpp_filt()'s objective is select the most preferable given a list of bn
#'   where bn list is specified in 2xN matrix that the earlier column the more preferable
#'
#'  gbp2d_solver_dpp_filt() use an approx binary search and determine f w.r.t bn.n_cols
#'   where f = 1 indicate the bn being selected and only one of 1 in result returned.
#'
#'  ok = true if any bin can fit all it and algorithm will select smallest bn can fit all
#'   otherwise ok = false and algorithm will select a bn can maximize volume of fitted it
#'
#'  often recommend to make the last and least preferable bn dominate all other bn in list
#'   when design bn list, bnX dominant bnY if all(X(l, d) > Y(l, d)).
#' @param ld
#'  it scale <matrix>
#'  - l, d it scale along x and y <numeric>
#' @param m
#'  bn scale <matrix>
#'  - l, d bn scale along x and y <numeric>
#'  - l, d in row and each col is a single bn
#'
#'  should make sure bn list are sorted via volume
#'   so that the first col is the most prefered smallest bn, and also
#'   the last col is the least prefered largest and often dominant bn
#'
#'  should make sure no X in front of Y if bnX dominant bnY,
#'   bnX dominant bnY if all(X(l, d) > Y(l, d)) and should always prefer Y.
#'
#'  should make sure bn such that l >= d or vice versa.
#' @return gbp2q
#'  a gbp2q instantiate with p profit, it item (x, y, l, d) position scale matrix, bn bin (l, d) scale matrix,
#'   k it selection, o objective, f bn selection, and ok an indicator of all fit or not.
#' @family gbp2q
#' @rdname gbp2d_solver_dpp_filt
"gbp2d_solver_dpp_filt"


#' gbp2q_checkr
#' @description
#'  auxilium of gbp2q and gbp2d_solver_dpp_filt
#' @details
#'  check fit solution is valid:
#'   no conflict between item and bin, and no conflict between each pair of item.
#' @param sn <gbp2q>
#'  gbp2q object from gbp2d_solver_dpp_filt() solution.
#' @return okfit? <bool>
#' @family gbp2q
#' @rdname gbp2q_checkr
"gbp2q_checkr"

