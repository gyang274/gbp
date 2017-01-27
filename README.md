<!-- README.md is generated from README.Rmd. Please edit that file -->
**gbp: A bin packing problem solver**
=====================================

Overview
--------

Basic infrastructure and several algorithms for **1d - 4d bin packing problem**. The **gbp** package provides a set of **c-level classes and solvers for 1d - 4d bin packing problem**, and **an r-level solver for 4d bin packing problem**, which is a wrapper over the c-level 4d bin packing problem solver.

The 4d bin packing problem solver aims to solve bin packing problem, a.k.a container loading problem, with an additional constraint on weight. Given a set of rectangular-shaped items, and a set of rectanular-shaped bins with weight limit, the solver looks for an orthogonal packing solution such that minimizes the number of bins and maximize volume utilization. Each rectangular-shaped item i = 1, .. , n is characterized by length l\_i, depth d\_i, height h\_i, and weight w\_i, and each rectanular-shaped bin j = 1, .. , m is specified similarly by length l\_j, depth d\_j, height h\_j, and weight limit w\_j. The item can be rotated into any orthogonal direction, and no further restrictions implied.

Vignettes
---------

TODO: add link to vignettes

[Shiny App](https://gyang.shinyapps.io/gbp_app/)
------------------------------------------------

A [shiny application](https://gyang.shinyapps.io/gbp_app/) that demonstrate how to use how to use **gbp** function **bpp\_solver** in fulfilling the order packing process in business operations.

OpenCPU API
-----------

An example call to **gbp::bpp\_solver** when an OpenCPU server is running with package **gbp** installed.

    curl http://localhost:8004/ocpu/library/gbp/R/bpp_solver/json -H "Content-Type: application/json" -d '{"it":[{"oid":1000001,"sku":"A0A0A0","l":2.14,"d":3.58,"h":4.76,"w":243},{"oid":1000001,"sku":"A0A0A1","l":7.24,"d":7.24,"h":2.58,"w":110},{"oid":1000001,"sku":"A0A0A1","l":7.24,"d":7.24,"h":2.58,"w":110},{"oid":1000002,"sku":"A0A0A0","l":2.14,"d":3.58,"h":4.76,"w":243},{"oid":1000002,"sku":"A0A0A1","l":7.24,"d":7.24,"h":2.58,"w":110},{"oid":1000002,"sku":"A0A0A1","l":7.24,"d":7.24,"h":2.58,"w":110},{"oid":1000002,"sku":"A0A0A2","l":6,"d":6,"h":6,"w":235},{"oid":1000002,"sku":"A0A0A3","l":4,"d":4,"h":4,"w":258}], "bn":[{"id":"K0001","l":10,"d":10,"h":10},{"id":"K0010","l":20,"d":20,"h":20}], "wlmt":800, "wprecision":1}'
