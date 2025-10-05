The `dune-foamgrid` module
==========================

The `dune-foamgrid` module is an implementation of the `dune-grid` interface that implements one- and two-dimensional grids in a physical space of arbitrary dimension. The grids are not expected to have a manifold structure, i.e., more than two elements can share a common facet. This makes FoamGrid the grid data structure of choice for simulating structures such as foams, discrete fracture networks, or network flow problems.

Element Parametrizations
------------------------

`dune-foamgrid` allows to use element parametrizations to improve the geometry approximation for curved domains. Each coarse grid element can be given a parametrization that describes an embedding into physical space. This does not influence the grid itself -- elements of a FoamGrid are always affine. However, when refining the grid, the new vertices are determined according to the parametrization. That way, the grid approaches the shape described by the parametrization functions more and more as it gets refined.

Grid Growth
-----------

As a unique feature of a `dune-grid` implementation, a FoamGrid is allowed to grow and shrink, i.e., elements can be added to or removed from the grid at runtime. Data attached to the grid is not invalidated by this. This allows to simulate problems with growing and shrinking domains, which is useful for various network growth and remodeling problems.

Installation
------------

`dune-foamgrid` requires the DUNE core modules, version 2.4 or later.

Please see the [general instructions for building DUNE modules](https://www.dune-project.org/doc/installation-notes.html) for detailed instructions on how to build the module.

Development
-----------

The [development version of `dune-foamgrid`](https://gitlab.dune-project.org/extensions/dune-foamgrid) can be obtained from the DUNE project's Gitlab installation.
At the same place an [issue tracker](https://gitlab.dune-project.org/extensions/dune-foamgrid/issues) can be found.

Publications
------------

The concepts of `dune-foamgrid` are presented in the following publication:

[O. Sander, T. Koch, N. Schröder, B. Flemisch, The Dune FoamGrid  implementation for surface and network grids, arXiv:1511.03415, 2015](http://arxiv.org/abs/1511.03415)
