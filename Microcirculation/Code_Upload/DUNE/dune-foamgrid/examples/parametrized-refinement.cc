// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/** \file
 * \brief Construct a FoamGrid with element parametrizations given by a closed-form function
 */

#include <config.h>
#include <iostream>
#include <cmath>
#include <memory>

#include <dune/common/function.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/foamgrid/foamgrid.hh>

using namespace Dune;

/**
 * \brief Mapping from R^d to the graph of a given function
 */
template<int dim, int dimworld>
class GraphMapping :
  public Dune::VirtualFunction<Dune::FieldVector<double, dim>, Dune::FieldVector<double, dimworld> >
{
  // Element corners in the global parameter domain
  std::array<FieldVector<double,dim>, dim+1> corners_;

  std::function<FieldVector<double,3>(FieldVector<double,2>)> graph_;

public:
  GraphMapping(std::array<FieldVector<double,dim>, dim+1> corners,
               std::function<FieldVector<double,3>(FieldVector<double,2>)> graph)
  : corners_(corners), graph_(graph)
  {}

  /**
   * \brief Function evaluation.
   *
   * \param x Argument for function evaluation.
   * \param y Result of function evaluation.
   */
  void evaluate(const Dune::FieldVector<double,dim>& x, Dune::FieldVector<double,dimworld>& y) const
  {
    // Linear interpolation between the corners
    auto globalX = corners_[0];
    for (size_t i=0; i<x.size(); i++)
      for (int j=0; j<dim; j++)
        globalX[j] += x[i]*(corners_[i+1][j]-corners_[0][j]);
    y = graph_(globalX);
  }

};

int main (int argc, char *argv[]) try
{
  const int dim      = 2;
  const int dimworld = 3;
  bool adaptive      = false;

  // Global parametrization function
  auto parametrization = [](const FieldVector<double,2>& x) -> FieldVector<double,3>
                           {return {x[0], x[1], 0.2*exp(-x.two_norm())*cos(4.5*M_PI*x.two_norm())};};


  // Create the grid
  typedef Dune::FoamGrid<dim, dimworld> Grid;

  // Start grid creation
  Dune::GridFactory<Grid> factory;

  // The list of grid vertex positions
  std::vector<FieldVector<double,2> > vertices = {{-1,-1}, {1,-1}, {-1, 1}, {1, 1}};

  for (const auto& p : vertices)
    factory.insertVertex(parametrization(p));

  // Create the element geometries
  Dune::GeometryType triangle;
  triangle.makeTriangle();

  std::array<FieldVector<double,2>, 3> corners0 = {vertices[0], vertices[1], vertices[2]};
  std::array<FieldVector<double,2>, 3> corners1 = {vertices[1], vertices[3], vertices[2]};

  factory.insertElement(triangle, {0,1,2}, std::make_shared<GraphMapping<dim, dimworld> >(corners0, parametrization));
  factory.insertElement(triangle, {1,3,2}, std::make_shared<GraphMapping<dim, dimworld> >(corners1, parametrization));

  // create the grid
  auto grid = factory.createGrid();

  // output VTK
  Dune::VTKWriter<Grid::LeafGridView > writer(grid->leafGridView());

  for (int i=0; i<6; i++)
  {
    writer.write("refine-" + std::to_string(i));
    for (const auto& element : elements(grid->leafGridView()))
    {
      auto center = element.geometry().center();
      if (adaptive == false or i==0 or std::fabs(center[1]) < std::pow(0.5,i-1) or center[1]>0)
        grid->mark(1,element);
    }

    grid->preAdapt();
    grid->adapt();
    grid->postAdapt();
  }
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Dune::Exception e) {
  std::cout << e << std::endl;
  return 1;
}
