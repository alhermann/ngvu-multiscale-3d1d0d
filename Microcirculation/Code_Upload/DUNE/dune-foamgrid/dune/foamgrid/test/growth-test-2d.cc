// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#include <config.h>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/test/gridcheck.hh>
#include <dune/grid/test/checkadaptation.hh>
#include <dune/grid/test/checkindexset.hh>

template <class Grid>
void checkGridElementGrowth(Grid& grid)
{
  using namespace Dune;
  enum { dimworld = Grid::dimensionworld };
  enum { dim = Grid::dimension };

  // the vertex mapper
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<Grid,Dune::MCMGVertexLayout> mapper(grid);

  // let the top element grow
  for (const auto& element : elements(grid.leafGridView()))
    for (const auto& intersection : intersections(grid.leafGridView(), element))
      if(intersection.boundary() && intersection.geometry().center()[1] > 0.5)
      {
        const Dune::ReferenceElement<double,dim>& refElement
                    = Dune::ReferenceElements<double,dim>::general(element.type());
        auto&& v0 = element.template subEntity<dim>(refElement.subEntity(intersection.indexInInside(), dim-1, 0, dim));
        auto&& v1 = element.template subEntity<dim>(refElement.subEntity(intersection.indexInInside(), dim-1, 1, dim));

        // calculate new vertex position
        Dune::FieldVector<double, dimworld> newVertex = v0.geometry().center();
        newVertex += v1.geometry().center(); newVertex /= 2.0;
        newVertex[1] += 0.5;

        // insert new vertex
        auto vIdx = grid.insertVertex(newVertex);

        // insert the element
        grid.insertElement(element.type(), {mapper.index(v0), mapper.index(v1), vIdx});
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
    DUNE_THROW(InvalidStateException,"grid.grow() does not return correct information");

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
        std::cout << "Boundary Intersection no"<<isCounter<<" has segment index: " << intersection.boundarySegmentIndex() << std::endl;
        ++isCounter;
      }
    }
  }

  numBoundarySegments = grid.numBoundarySegments();
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << "numBoundarySegments after growth: " << numBoundarySegments << std::endl<< std::endl;
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

using namespace Dune;

int main (int argc, char *argv[])
{
  try
  {
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    static const int dimworld = 3;
    static const int dim = 2;
    std::cout << "Creating a FoamGrid<"<<dim<<", "<<dimworld<<"> ("<<dim<<"d in "<<dimworld<<"d grid)" << std::endl;
    Dune::FieldVector<double,dimworld> lower(0.0);
    Dune::FieldVector<double,dimworld> upper(1.0);
    std::array<unsigned int,dim> elements;
    std::fill(elements.begin(), elements.end(), 1);
    std::shared_ptr<FoamGrid<dim,dimworld> > grid = StructuredGridFactory<FoamGrid<dim,dimworld> >::createSimplexGrid(lower,upper,elements);

    Dune::VTKWriter<typename FoamGrid<dim,dimworld>::LeafGridView > writer(grid->leafGridView(), VTK::nonconforming);
    writer.write("before_growth");

    // check simple grid growth
    Dune::gridinfo(*grid);
    checkGridElementGrowth(*grid);
    writer.write("after_growth");
    Dune::gridinfo(*grid);

    // check removal of a grid element
    checkGridElementRemoval(*grid);
    writer.write("after_removal");
    Dune::gridinfo(*grid);
    checkIndexSet(*grid, grid->leafGridView(), std::cout);
    for (int i = 0; i < grid->maxLevel(); ++i)
        checkIndexSet(*grid, grid->levelGridView(i), std::cout);

    // add an element unrelated to the rest of the grid
    auto idx1 = grid->insertVertex({3.0, 3.0, 3.0});
    auto idx2 = grid->insertVertex({3.0, 3.5, 3.0});
    auto idx3 = grid->insertVertex({3.5, 3.0, 3.0});
    grid->insertElement(Dune::GeometryType(Dune::GeometryType::simplex, dim), {idx1, idx2, idx3});
    grid->preGrow();
    grid->grow();
    grid->postGrow();
    writer.write("add_unrelated_element");
    Dune::gridinfo(*grid);
    checkIndexSet(*grid, grid->leafGridView(), std::cout);
    for (int i = 0; i < grid->maxLevel(); ++i)
      checkIndexSet(*grid, grid->levelGridView(i), std::cout);

    // call growth without adding anything
    bool removedElements = grid->preGrow();
    if(removedElements)
      DUNE_THROW(InvalidStateException,"grid.preGrow() does not return correct information");
    bool addedElements = grid->grow();
    if(addedElements)
      DUNE_THROW(InvalidStateException,"grid.grow() does not return correct information");
    grid->postGrow();

    // do a grid check on a refined grid
    grid->globalRefine(4);
    gridcheck(*grid);
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
