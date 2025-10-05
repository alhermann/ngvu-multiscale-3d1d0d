#install headers
install(FILES foamgrid.cc
              foamgridedge.hh
              foamgridelements.hh
              foamgridentity.hh
              foamgridentitypointer.hh
              foamgridentityseed.hh
              foamgridfactory.hh
              foamgridgeometry.hh
              foamgridhierarchiciterator.hh
              foamgridindexsets.hh
              foamgridintersectioniterators.hh
              foamgridintersections.hh
              foamgridleafiterator.hh
              foamgridleveliterator.hh
              foamgridnulliteratorfactory.hh
              foamgridvertex.hh
              foamgridviews.hh
         DESTINATION include/dune/foamgrid/foamgrid)
