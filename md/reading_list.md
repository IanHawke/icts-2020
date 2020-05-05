# Reading list

This is a very incomplete list of sources I regularly use.

#### Reviews

* [Font, Numerical Hydrodynamics and Magnetohydrodynamics in General Relativity, Living Review](https://doi.org/10.12942/lrr-2008-7). Last updated in 2008 but crucial background.
* [Marti & Müller, Grid-based Methods in Relativistic Hydrodynamics and Magnetohydrodynamics, Living Review](https://doi.org/10.1007/lrca-2015-3). SR only but updated in 2015. See also the [original version](https://link.springer.com/article/10.12942/lrr-2003-7) from 2003.
* [Balsara, Higher-order accurate space-time schemes for computational astrophysics—Part I: finite volume methods, Living Review](https://doi.org/10.1007/s41115-017-0002-8). Very methods heavy. Cutting edge but not easy going.
* [Baiotti & Rezzolla, Binary neutron star mergers: a review of Einstein's richest laboratory](https://doi.org/10.1088/1361-6633/aa67bb) ([arxiv:1607.03540](https://arxiv.org/pdf/1607.03540.pdf)). From 2017, specific to mergers.
* [Shibata & Taniguchi, Coalescence of Black Hole-Neutron Star Binaries](https://link.springer.com/article/10.12942/lrr-2011-6). From 2011, so somewhat dated: follow up by looking at Shibata's lengthy track-record in binary merger simulations.

#### Theses

* [Radice, Advanced Numerical Approaches in the Dynamics of Relativistic Flows](https://www.astro.princeton.edu/~dradice/downloads/thesis.pdf). From 2013, touches on a number of important technical details.

#### Books

* [Leveque, Finite Volume Methods for Hyperbolic Problems, CUP](https://www.cambridge.org/core/books/finite-volume-methods-for-hyperbolic-problems/97D5D1ACB1926DA1D4D52EAD6909E2B9). No astrophysics but one of the standard numerical methods texts.
* [Hesthaven, Numerical Methods for Conservation Laws: From Analysis to Algorithms, SIAM](https://epubs.siam.org/doi/book/10.1137/1.9781611975109). Still no astrophysics and even more mathematical-technical, but goes deep into methods like Discontinuous Galerkin and spectral elements which may be the future direction of the field.
* [Rezzolla & Zanotti, Relativistic Hydrodynamics, OUP](https://global.oup.com/academic/product/relativistic-hydrodynamics-9780198528906). From 2013, its focus is on hydrodynamics, not MHD. Lots of detail.
* [Alcubierre, Introduction to 3+1 Numerical Relativity, OUP](https://global.oup.com/academic/product/introduction-to-31-numerical-relativity-9780199656158). From 2012, its focus is really vacuum relativity, but introduces hydrodynamics well from that viewpoint.
* [Andersson, Gravitational-Wave Astronomy, OUP](https://global.oup.com/academic/product/gravitational-wave-astronomy-9780198568032). Much more on the neutron star modelling, with some chapters on where the numerics fits in. From 2019.

#### Codes and tutorials

* [Open Astrophysics Bookshelf](https://open-astrophysics-bookshelf.github.io/). Relativity isn't a focus but the material covers a lot of numerics in great depth, with example codes throughout. Have a look at [github.com/python-hydro](https://github.com/python-hydro) for detailed examples in one and two dimensions.
* [NRPy](https://github.com/zachetienne/nrpytutorial). A Python front end to a numerical relativity code, linked to [BlackHoles@Home](http://astro.phys.wvu.edu/bhathome/). Focus is on vacuum, but there are examples from GR(M)HD.
* [Einstein Toolkit](https://einsteintoolkit.org/). A production GRMHD code that runs on massively parallel machines. There's a steep learning curve and it's designed to do a broad range of things (so it's more complex than it needs to be to do any *one* thing), but this can be used for real research.
