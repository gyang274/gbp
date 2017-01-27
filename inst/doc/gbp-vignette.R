## ----message=FALSE-------------------------------------------------------
library(gbp)

## ------------------------------------------------------------------------
#- bpp_solver: input: order list in data.table `it` and bin list in data.table `bn`

#- it
#  it item <data.table>
#  - oid order id <integer>
#  - sku items id <character>
#  - l it length which scale will be placed along x-coordinate <numeric>
#  - d it depth  which scale will be placed along y-coordinate <numeric>
#  - h it height which scale will be placed along z-coordinate <numeric>
#  - w it weight which scale will be placed along w-coordinate <numeric>
# l d h are subject to rotate, while w is on a separate single dimension
it <- data.table::data.table(
  oid = c(1428571L, 1428571L, 1428571L, 1428572L, 1428572L, 1428572L, 1428572L, 1428572L),
  sku = c("A0A0A0", "A0A0A1", "A0A0A1", "A0A0A0", "A0A0A1", "A0A0A1", "A0A0A2", "A0A0A3"),
  l   = c(2.140000, 7.240000, 7.240000, 2.140000, 7.240000, 7.240000, 6.000000, 4.000000),
  d   = c(3.580000, 7.240000, 7.240000, 3.580000, 7.240000, 7.240000, 6.000000, 4.000000),
  h   = c(4.760000, 2.580000, 2.580000, 4.760000, 2.580000, 2.580000, 6.000000, 4.000000),
  w   = c(243.0000, 110.0000, 110.0000, 243.0000, 110.0000, 110.0000, 235.0000, 258.0000)
)

knitr::kable(it)

#- bn
#  bn bins <data.table>
#  - id bn id <character>
#  - l: bn length limit along x-coordinate <numeric>
#  - d: bn depth  limit along y-coordinate <numeric>
#  - h: bn height limit along z-coordinate <numeric>
#  - w: bn weight limit along w - a separate single dimension <numeric>
#  - l, d, h will be sorted to have l >= d >= h within solver
# bin must be ordered by preference such that the first bin is most preferred one.

bn <- data.table::data.table(
  id = c("K0001", "K0002", "K0003", "K0004", "K0005"),
  l  = c(06.0000, 10.0000, 09.0000, 10.0000, 22.0000),
  d  = c(06.0000, 08.0000, 08.0000, 10.0000, 14.0000),
  h  = c(06.0000, 06.0000, 07.0000, 10.0000, 09.0000),
  w  = c(600.000, 600.000, 800.000, 800.000, 800.000)
)


knitr::kable(bn)


## ------------------------------------------------------------------------
#- bpp_solver: output: packing solution

#- sn
#  sn bpp_solution packing solution <list>
#  - it item <data.table>
#    - oid: order id <integer>
#    - sku: stock keeping unit as it id <character>
#    - tid: ticket id - an unique id within oid <integer>
#    - otid: order id x ticket id - an unique indentifier indicate it with same tid can be packed into one bin <character>
#    - bid: bn id <integer>
#    - x, y, z it position in bid bin <numeric>
#    - l, d, h it scale along x, y, z <numeric>
#    - w it weight <numeric>
#  - bn bins <data.table>
#    - id bn id <character>
#    - l bn length limit along x-coordinate <numeric>
#    - d bn depth  limit along y-coordinate <numeric>
#    - h bn height limit along z-coordinate <numeric>
#    - w bn weight limit along w - a separate single dimension <numeric>
sn <- gbp::bpp_solver(it = it, bn = bn)

sn$it

sn$bn


## ------------------------------------------------------------------------
#- bpp_viewer: a packing solution viewer
#   requires rgl, run in r-console will create interactive 3d views in rgl window
# bpp_viewer(sn)

## ------------------------------------------------------------------------
#- gbp4d

#- ldhw: item l, d, h, w in matrix
ldhw <- t(as.matrix(it[oid == 1428572L, .(l, d, h, w)]))
ldhw

#- m: bin l, d, h in matrix
m <- t(as.matrix(bn[ , .(l, d, h, w)])) # multple bin
m

#- p: item fit sequence w.r.t bin
p <- gbp4d_solver_dpp_prep_create_p(ldhw, m[, 4L]) # single bin
p

#- sn
sn4d <- gbp4d_solver_dpp(p, ldhw, m[, 4L])

sn4d$it # matrix of items x, y, z, w (weight bn is holding when fit it into bn), l, d, h, w (weight of it itself) (x, y, z, w set into -1 when item not fit into bin)

sn4d$k  # indicator of which items are fitted into bin

sn4d$bn # matrix of bins l, d, h, w (weight limit)

sn4d$o  # volume of items fitted into bin

sn4d$ok # indicator of whether all items are fitted into bin

## ------------------------------------------------------------------------
gbp4d_checkr(sn4d) #- check: no overlap in 3d volume and no over weight limit

## ------------------------------------------------------------------------
# gbp4d_viewer(sn4d)

## ------------------------------------------------------------------------
#- gbp4q

sn4q <- gbp4d_solver_dpp_filt(ldhw, m) # multiple bins, no fit sequence p 
# p is determined w.r.t each bin using algorithm in gbp4d_solver_dpp_prep_create_p 

sn4q$it # matrix of items x, y, z, w (weight bn is holding when fit it into bn), l, d, h, w (weight of it itself) (x, y, z, w set into -1 when item not fit into bin)

sn4q$k  # indicator of which items are fitted into bin

sn4q$bn # matrix of bins l, d, h, w (weight limit)

sn4q$f  # indicator of which bin is selected

sn4q$o  # volume of items fitted into bin

sn4q$ok # indicator of whether all items are fitted into bin

## ------------------------------------------------------------------------
gbp4q_checkr(sn4q) #- check: no overlap in 3d volume and no over weight limit

## ------------------------------------------------------------------------
# gbp4q_viewer(sn4q)

## ------------------------------------------------------------------------
#- gbp3d

sn3d <- gbp3d_solver_dpp(p, ldhw[1L:3L, ], m[, 4L])

sn3d$it # matrix of items x, y, z, l, d, h (x, y, z set into -1 when item not fit into bin)

sn3d$k  # indicator of which items are fitted into bin

sn3d$bn # matrix of bins l, d, h

sn3d$o  # volume of items fitted into bin

sn3d$ok # indicator of whether all items are fitted into bin

## ------------------------------------------------------------------------
gbp3d_checkr(sn3d) #- check: no overlap in 3d volume

## ------------------------------------------------------------------------
# gbp3d_viewer(sn3d)

## ------------------------------------------------------------------------
#- gbp3q

sn3q <- gbp3d_solver_dpp_filt(ldhw[1L:3L, ], m) # multiple bins, no fit sequence p
# p is determined w.r.t each bin using algorithm in gbp3d_solver_dpp_prep_create_p 

sn3q$it # matrix of items x, y, z, l, d, h (x, y, z set into -1 when item not fit into bin)

sn3q$k  # indicator of which items are fitted into bin

sn3q$bn # matrix of bins l, d, h

sn3q$f  # indicator of which bin is selected

sn3q$o  # volume of items fitted into bin

sn3q$ok # indicator of whether all items are fitted into bin

## ------------------------------------------------------------------------
gbp3q_checkr(sn3q) #- check: no overlap in 3d volume

## ------------------------------------------------------------------------
# gbp3q_viewer(sn3q)

## ------------------------------------------------------------------------
#- gbp2d

sn2d <- gbp2d_solver_dpp(p, ldhw[1L:2L, ], m[, 4L])

sn2d$it # matrix of items x, y, l, d (x, y set into -1 when item not fit into bin)

sn2d$k  # indicator of which items are fitted into bin

sn2d$bn # matrix of bins l, d

sn2d$o  # volume of items fitted into bin

sn2d$ok # indicator of whether all items are fitted into bin

## ------------------------------------------------------------------------
gbp2d_checkr(sn2d) #- check: no overlap in 2d area

## ------------------------------------------------------------------------
# gbp2d_viewer(sn2d) #- view on XZ surface and set Y into 1 to give stereo perception

## ------------------------------------------------------------------------
#- gbp2q

sn2q <- gbp2d_solver_dpp_filt(ldhw[1L:2L, ], m) # multiple bins, no fit sequence p
# p is determined w.r.t each bin using algorithm in gbp2d_solver_dpp_prep_create_p 

sn2q$it # matrix of items x, y, l, d (x, y set into -1 when item not fit into bin)

sn2q$k  # indicator of which items are fitted into bin

sn2q$bn # matrix of bins l, d

sn2q$f  # indicator of which bin is selected

sn2q$o  # volume of items fitted into bin

sn2q$ok # indicator of whether all items are fitted into bin

## ------------------------------------------------------------------------
gbp2q_checkr(sn2q) #- check: no overlap in 2d area

## ------------------------------------------------------------------------
# gbp2q_viewer(sn2q)

## ------------------------------------------------------------------------
#- gbp1d

v <- apply(ldhw[1L:3L, ], 2L, prod)

sn1d <- gbp1d_solver_dpp(p = v, w = ldhw[4L, ], c = 714.28)

sn1d$p # vector of items profit

sn1d$w # vector of items weight

sn1d$c # weight limit constraint

sn1d$k  # indicator of which items are selected

sn1d$o  # weight of items selected

sn1d$ok # indicator of whether all items are selected

