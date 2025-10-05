// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#include <config.h>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/test/gridcheck.hh>
#include <dune/grid/test/checkadaptation.hh>
#include <dune/grid/test/checkindexset.hh>

template <class Grid>
void checkGridElementGrowth(Grid& grid, int numElements)
{
  using namespace Dune;
  enum { dimworld = Grid::dimensionworld };
  enum { dim = Grid::dimension };

  for (const auto& element : elements(grid.leafGridView()))
  {
    const auto geo = element.geometry();
    if(geo.center()[1] >= 0.0 && geo.center()[0] < -1.0)
    {
      Dune::FieldVector<double, dimworld> growPoint = geo.center();
      double interval = 2.0/numElements;
      unsigned int oldVIdx = 0;
      for (int i=0; i<numElements; ++i)
      {
        growPoint += Dune::FieldVector<double, dimworld>(interval);

        // insert a new vertex
        auto newVIdx = grid.insertVertex(growPoint);

        if (i==0)
        {
          // find the second index
          Dune::LeafMultipleCodimMultipleGeomTypeMapper<Grid,Dune::MCMGVertexLayout> mapper(grid);

          // insert the new element
          grid.insertElement(Dune::GeometryType(1),
                             {mapper.index(element.template subEntity<dim>(0)), newVIdx}
                             );
        }
        else
        {
          grid.insertElement(Dune::GeometryType(1),
                             {oldVIdx, newVIdx}
                             );
        }
        oldVIdx = newVIdx;
      }
    }
  }

  std::size_t numBoundarySegments = grid.numBoundarySegments();
  std::cout << std::endl<< "numBoundarySegments before growth: " << numBoundarySegments << std::endl;
  std::cout << "-------------------------------------------" << std::endl;

  bool elementsWillVanish = grid.preGrow();
  if(elementsWillVanish)
    DUNE_THROW(InvalidStateException,"grid.preGrow() does not return correct information");

  // check mightVanish
  for (const auto& element : elements(grid.levelGridView(0)))
    checkHierarchy(element);

  bool newElementGenerated = grid.grow();
  if(!newElementGenerated)
    DUNE_THROW(InvalidStateException,"grid.preGrow() does not return correct information");

  // obtain growth insertion index of new elements
  for (const auto& element : elements(grid.levelGridView(0)))
  {
    if (element.isNew())
      std::cout << "After calling grow(): found new element with growth insertion index: " << grid.growthInsertionIndex(element) << std::endl;
  }

  grid.postGrow();

  // Loop over all levels except the lowest one
  for (int level = 0 ; level <= grid.maxLevel(); ++level )
  {
      for (const auto& element : elements(grid.levelGridView(level)))
        if(element.isNew())
        DUNE_THROW(InvalidStateException,"After postGrow() was called no entity is new, i.e., isNew() == false");
  }

  std::cout << "Boundary intersections after growth: " << std::endl;
  int isCounter = 0;
  for (const auto& element : elements(grid.leafGridView()))
  {
    for (const auto& intersection : intersections(grid.leafGridView(), element))
    {
      if(intersection.boundary())
      {
        std::cout << "--Boundary Intersection no"<<isCounter<<" has segment index: " << intersection.boundarySegmentIndex()
                  << " --> at position: " << intersection.geometry().center() << std::endl;
        ++isCounter;
      }
    }
  }

  numBoundarySegments = grid.numBoundarySegments();
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << "numBoundarySegments after growth: " << numBoundarySegments << std::endl<< std::endl;
}

template <class Grid>
void checkGridElementMerge(Grid& grid)
{
  enum { dim = Grid::dimension };

  // vertex mapper
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<Grid,Dune::MCMGVertexLayout> mapper(grid);

  // insert an element between the two closest boundary facets
  std::vector<unsigned int> vertices(dim+1);
  double dist = std::numeric_limits<double>::max();
  for (const auto& element : elements(grid.leafGridView()))
  {
    for (const auto& intersection : intersections(grid.leafGridView(), element))
    {
      auto center = intersection.geometry().center();
      if(!intersection.boundary())
        continue;
      for (const auto& element2 : elements(grid.leafGridView()))
      {
        if(element2 == element)
          continue;

        for (const auto& intersection2 : intersections(grid.leafGridView(), element2))
        {
          if(!intersection2.boundary())
            continue;
          auto center2 = intersection2.geometry().center();
          auto diff = center;
          diff -= center2;
          if(diff.two_norm() < dist)
          {
            dist = diff.two_norm();
            vertices[0] = mapper.index(element.template subEntity<dim>(intersection.indexInInside()));
            vertices[1] = mapper.index(element2.template subEntity<dim>(intersection2.indexInInside()));
          }
        }
      }
    }
  }
  grid.insertElement(Dune::GeometryType(1), vertices);

  std::size_t numBoundarySegments = grid.numBoundarySegments();
  std::cout << std::endl<< "numBoundarySegments before merge: " << numBoundarySegments << std::endl;
  std::cout << "-------------------------------------------" << std::endl;

  bool newElementInserted = grid.grow();
  if(!newElementInserted)
    DUNE_THROW(Dune::InvalidStateException,"grid.grow() does not return correct information");
  grid.postGrow();

  std::cout << "Boundary intersections after merge: " << std::endl;

  int isCounter = 0;
  for (const auto& element : elements(grid.leafGridView()))
  {
    for (const auto& intersection : intersections(grid.leafGridView(), element))
    {
      if(intersection.boundary())
      {
        std::cout << "Boundary Intersection no"<<isCounter<<" has segment index: " << intersection.boundarySegmentIndex() << std::endl;
        ++isCounter;
      }
    }
  }

  numBoundarySegments = grid.numBoundarySegments();
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << "numBoundarySegments after merge: " << numBoundarySegments << std::endl<< std::endl;
}

template <class Grid>
void checkGridElementRemoval(Grid& grid)
{
  using namespace Dune;
  enum { dimworld = Grid::dimensionworld };
  enum { dim = Grid::dimension };

  // remove the first inserted element
  auto eIt = grid.leafGridView().template begin<0>();
  grid.removeElement(*eIt);

  bool elementsWillVanish = grid.preGrow();
  if(!elementsWillVanish)
    DUNE_THROW(InvalidStateException,"grid.preGrow() does not return correct information");

  std::size_t numBoundarySegments = grid.numBoundarySegments();
  std::cout << std::endl << "numBoundarySegments before removal: " << numBoundarySegments << std::endl;
  std::cout << "-------------------------------------------" << std::endl;

  bool newElementInserted = grid.grow();
  if(newElementInserted)
    DUNE_THROW(Dune::InvalidStateException,"grid.grow() does not return correct information");
  grid.postGrow();

  std::cout << "Boundary intersections after removal: " << std::endl;

  int isCounter = 0;
  for (const auto& element : elements(grid.leafGridView()))
  {
    for (const auto& intersection : intersections(grid.leafGridView(), element))
    {
      if(intersection.boundary())
      {
        std::cout << "Boundary Intersection no"<<isCounter<<" has segment index: " << intersection.boundarySegmentIndex() << std::endl;
        ++isCounter;
      }
    }
  }

  numBoundarySegments = grid.numBoundarySegments();
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << "numBoundarySegments after removal: " << numBoundarySegments << std::endl << std::endl;
}

template <class Grid>
void checkGridElementGrowthLevel(Grid& grid)
{
  enum { dim = Grid::dimension };
  enum { dimworld = Grid::dimensionworld };

  // first we partially refine the grid at the xMax boundary
  double xMax = std::numeric_limits<double>::min();
  double xMin = std::numeric_limits<double>::max();
  for (const auto& vertex : vertices(grid.leafGridView()))
  {
    xMax = std::max(xMax, vertex.geometry().center()[0]);
    xMin = std::min(xMin, vertex.geometry().center()[0]);
  }

  for (const auto& element : elements(grid.leafGridView()))
    for (const auto& intersection : intersections(grid.leafGridView(), element))
      if(intersection.boundary())
        if(intersection.geometry().center()[0] == xMax)
          grid.mark(1, element);

  grid.preAdapt();
  grid.adapt();
  grid.postAdapt();

  std::cout << std::endl << "-------------------------------------------" << std::endl;
  std::cout << "gridinfo after adaptation: " << std::endl;
  Dune::gridinfo(grid);
  std::cout << "-------------------------------------------" << std::endl;

  // vertex mapper
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<Grid,Dune::MCMGVertexLayout> mapper(grid);
  // Then we do a merge with a level 1 and a level 0 vertex
  std::vector<unsigned int> elementVertices(dim+1);
  for (const auto& element : elements(grid.leafGridView()))
    for (const auto& intersection : intersections(grid.leafGridView(), element))
      if(intersection.boundary())
      {
        if(intersection.geometry().center()[0] == xMax)
        {
          auto vertex = element.template subEntity<dim>(intersection.indexInInside());
          assert(vertex.level() == 1);
          elementVertices[0] = mapper.index(vertex);
        }
        else if(intersection.geometry().center()[0] == xMin)
        {
          auto vertex = element.template subEntity<dim>(intersection.indexInInside());
          assert(vertex.level() == 0);
          elementVertices[1] = mapper.index(vertex);
        }
      }

  grid.insertElement(Dune::GeometryType(1), elementVertices);

  bool elementsWillVanish = grid.preGrow();
  if(elementsWillVanish)
    DUNE_THROW(Dune::InvalidStateException,"grid.preGrow() does not return correct information");

  std::size_t numBoundarySegments = grid.numBoundarySegments();
  std::cout << std::endl << "numBoundarySegments before complex merge: " << numBoundarySegments << std::endl;
  std::cout << "-------------------------------------------" << std::endl;

  bool newElementInserted = grid.grow();
  if(!newElementInserted)
    DUNE_THROW(Dune::InvalidStateException,"grid.grow() does not return correct information");

  for (const auto& element : elements(grid.leafGridView()))
    if(element.isNew())
      assert(element.level() == 0);

  grid.postGrow();

  std::cout << "Boundary intersections after complex merge: " << std::endl;

  int isCounter = 0;
  for (const auto& element : elements(grid.leafGridView()))
  {
    for (const auto& intersection : intersections(grid.leafGridView(), element))
    {
      if(intersection.boundary())
      {
        std::cout << "Boundary Intersection no"<<isCounter<<" has segment index: " << intersection.boundarySegmentIndex() << std::endl;
        ++isCounter;
      }
    }
  }

  numBoundarySegments = grid.numBoundarySegments();
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << "numBoundarySegments after complex merge: " << numBoundarySegments << std::endl << std::endl;

  Dune::FieldVector<double, dimworld> growPoint(6.0);
  growPoint[0] = 0.0;
  auto vIdx = grid.insertVertex(growPoint);
  grid.insertElement(Dune::GeometryType(1), {vIdx, elementVertices[0]});
  grid.insertElement(Dune::GeometryType(1), {vIdx, elementVertices[1]});

  elementsWillVanish = grid.preGrow();
  if(elementsWillVanish)
    DUNE_THROW(Dune::InvalidStateException,"grid.preGrow() does not return correct information");

  numBoundarySegments = grid.numBoundarySegments();
  std::cout << std::endl << "numBoundarySegments before complex merge: " << numBoundarySegments << std::endl;
  std::cout << "-------------------------------------------" << std::endl;

  newElementInserted = grid.grow();
  if(!newElementInserted)
    DUNE_THROW(Dune::InvalidStateException,"grid.grow() does not return correct information");

  for (const auto& element : elements(grid.leafGridView()))
    if(element.isNew())
      assert(element.level() == 0);

  grid.postGrow();

  std::cout << "Boundary intersections after complex merge: " << std::endl;

  isCounter = 0;
  for (const auto& element : elements(grid.leafGridView()))
  {
    for (const auto& intersection : intersections(grid.leafGridView(), element))
    {
      if(intersection.boundary())
      {
        std::cout << "Boundary Intersection no"<<isCounter<<" has segment index: " << intersection.boundarySegmentIndex() << std::endl;
        ++isCounter;
      }
    }
  }
}

using namespace Dune;

int main (int argc, char *argv[])
{
  try
  {
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    ////////////////////////////////1D-2D/////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    const std::string dune_foamgrid_path = std::string(DUNE_FOAMGRID_EXAMPLE_GRIDS_PATH) + "gmsh/";
    std::cout << "Creating a FoamGrid<1, 2> (1d in 2d grid)" << std::endl;
    std::shared_ptr<FoamGrid<1, 2> > grid1d( GmshReader<FoamGrid<1, 2> >::read(dune_foamgrid_path + "line1d2d.msh", /*verbose*/ true, false ) );

    Dune::VTKWriter<typename FoamGrid<1, 2>::LeafGridView > writer(grid1d->leafGridView(), VTK::nonconforming);
    writer.write("before_growth_1d2d");

    // check simple grid growth
    std::cout << std::endl << "Check simple grid growth (add element)" << std::endl;
    Dune::gridinfo(*grid1d);
    checkGridElementGrowth(*grid1d, 3);
    writer.write("after_growth_1d2d");
    Dune::gridinfo(*grid1d);

    // check a merger, i.e. inserting an element only with existing vertices
    std::cout << std::endl << "Check merger (connect two element facets)" << std::endl;
    checkGridElementMerge(*grid1d);
    writer.write("after_merge_1d2d");
    Dune::gridinfo(*grid1d);

    // check removal of a grid element
    std::cout << std::endl << "Check removal (remove element)" << std::endl;
    checkGridElementRemoval(*grid1d);
    writer.write("after_removal_1d2d");
    Dune::gridinfo(*grid1d);

    // check growth when vertices are on different levels
    std::cout << std::endl << "Check growth with vertices on different levels" << std::endl;
    checkGridElementGrowthLevel(*grid1d);
    writer.write("after_second_growth_1d2d");
    Dune::gridinfo(*grid1d);
    checkIndexSet(*grid1d, grid1d->leafGridView(), std::cout);
    for (int i = 0; i < grid1d->maxLevel(); ++i)
        checkIndexSet(*grid1d, grid1d->levelGridView(i), std::cout);

    // do a grid check on a refined grid
    std::cout << std::endl << "Check globalRefine on grown grid" << std::endl;
    grid1d->globalRefine(4);
    gridcheck(*grid1d);

    std::cout << std::endl << std::endl << std::endl;
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    ////////////////////////////////1D-3D/////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    std::cout << std::endl << std::endl << std::endl << std::endl;
    std::cout << "Creating a FoamGrid<1, 3> (1d in 3d grid)" << std::endl;
    std::shared_ptr<FoamGrid<1, 3> > grid1d3d( GmshReader<FoamGrid<1, 3> >::read(dune_foamgrid_path + "line1d2d.msh", /*verbose*/ true, false ) );

    Dune::VTKWriter<typename FoamGrid<1, 3>::LeafGridView > writer2(grid1d3d->leafGridView(), VTK::nonconforming);
    writer2.write("before_growth_1d3d");

    // check simple grid growth
    std::cout << std::endl << "Check simple grid growth (add element)" << std::endl;
    Dune::gridinfo(*grid1d3d);
    checkGridElementGrowth(*grid1d3d, 2);
    writer2.write("after_growth_1d3d");
    Dune::gridinfo(*grid1d3d);
    checkIndexSet(*grid1d3d, grid1d3d->leafGridView(), std::cout);

    // check a merger, i.e. inserting an element only with existing vertices
    std::cout << std::endl << "Check merger (connect two element facets)" << std::endl;
    checkGridElementMerge(*grid1d3d);
    writer2.write("after_merge_1d3d");
    Dune::gridinfo(*grid1d3d);

    // check removal of a grid element
    std::cout << std::endl << "Check removal (remove element)" << std::endl;
    checkGridElementRemoval(*grid1d3d);
    writer2.write("after_removal_1d3d");
    Dune::gridinfo(*grid1d3d);

    // check growth when vertices are on different levels
    std::cout << std::endl << "Check growth with vertices on different levels" << std::endl;
    checkGridElementGrowthLevel(*grid1d3d);
    writer2.write("after_second_growth_1d3d");
    Dune::gridinfo(*grid1d3d);
    checkIndexSet(*grid1d3d, grid1d3d->leafGridView(), std::cout);
    for (int i = 0; i < grid1d3d->maxLevel(); ++i)
        checkIndexSet(*grid1d3d, grid1d3d->levelGridView(i), std::cout);

    // do a grid check on a refined grid
    std::cout << std::endl << "Check globalRefine on grown grid" << std::endl;
    grid1d3d->globalRefine(4);
    gridcheck(*grid1d3d);

  }
  // //////////////////////////////////
  //   Error handler
  // /////////////////////////////////
  catch (Exception e) {
    std::cout << e << std::endl;
    return 1;
  }
  return 0;
}
