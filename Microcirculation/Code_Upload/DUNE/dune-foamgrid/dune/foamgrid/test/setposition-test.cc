#include <config.h>

#include <iostream>

#include <dune/grid/test/gridcheck.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/foamgrid/foamgrid.hh>


template<int dimworld>
Dune::FieldVector<double, dimworld>
getPosition(const Dune::FieldVector<double, dimworld>& pos)
{
  return Dune::FieldVector<double, dimworld>({pos[0], pos[1], pos[0]*pos[0] + pos[1]*pos[1]});
}

int main (int argc, char *argv[]) try
{
  // maybe initialize mpi
  Dune::MPIHelper::instance(argc, argv);

  std::cout << "Creating a FoamGrid<2, 3> (2d in 3d grid)" << std::endl;
  static const int dimworld = 3;
  static const int dimgrid = 2;
  typedef Dune::FoamGrid<dimgrid, dimworld> Grid;

  // build structured simplex grid
  Dune::FieldVector<double, dimworld> lower({-1, -1, 0});
  Dune::FieldVector<double, dimworld> upper({1, 1, 0});
  std::array<unsigned int, dimgrid> elements;
  std::fill(elements.begin(), elements.end(), 10);
  std::shared_ptr<Grid> grid =
    Dune::StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, elements);

  grid->globalRefine(2);

  Dune::VTKWriter<typename Grid::LeafGridView > writer(grid->leafGridView(), Dune::VTK::nonconforming);
  writer.write("setposition_0");

  // move all vertices
  for (const auto& vertex : vertices(grid->leafGridView()))
    grid->setPosition(vertex, getPosition(vertex.geometry().corner(0)));

  writer.write("setposition_1");

  grid->globalRefine(-1);

  writer.write("setposition_2");

  gridcheck(*grid);

  return 0;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Dune::Exception e) {
    std::cout << e << std::endl;
    return 1;
}
