// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#include <config.h>

#include "make2din3dgrid.hh"
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/test/gridcheck.hh>
#include <dune/grid/test/checkindexset.hh>
#include <dune/grid/test/checkintersectionit.hh>
#include <dune/grid/test/checkadaptation.hh>
#include <dune/grid/test/checkgeometryinfather.hh>
#include <dune/foamgrid/foamgrid.hh>

int main (int argc, char *argv[]) try
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    const std::string dune_foamgrid_path = std::string(DUNE_FOAMGRID_EXAMPLE_GRIDS_PATH) + "gmsh/";
    Dune::FieldVector<double,2> lower = {0,0};
    Dune::FieldVector<double,2> upper = {1,1};
    std::array<unsigned int,2> elements = {{6,6}};

    std::cout << "Checking FoamGrid<2, 2> (2d in 2d grid)" << std::endl;
    std::cout << "  Creating grid" << std::endl;
    std::shared_ptr<FoamGrid<2, 2> > grid2d = StructuredGridFactory<FoamGrid<2,2> >::createSimplexGrid(lower,upper,elements);
    {
        Dune::VTKWriter<typename FoamGrid<2, 2>::LeafGridView > writer(grid2d->leafGridView(), VTK::nonconforming);
        writer.write("2d_refined0");
    }
    gridcheck(*grid2d);
    Dune::gridinfo(*grid2d);

    // check grid adaptation interface
    checkAdaptation(*grid2d);

    // check grid after adaptive refinement
    checkAdaptRefinement(*grid2d);
    {
        Dune::VTKWriter<typename FoamGrid<2, 2>::LeafGridView > writer(grid2d->leafGridView(), VTK::nonconforming);
        writer.write("2d_refined-l");
    }
    Dune::gridinfo(*grid2d);
    checkGeometryInFather(*grid2d);
    gridcheck(*grid2d);
    checkIntersectionIterator(*grid2d);

    std::cout << "Checking FoamGrid<1, 3> (1d in 3d grid)" << std::endl;
    std::cout << "  Creating grid" << std::endl;
    std::shared_ptr<FoamGrid<1, 3> > grid1d( GmshReader<FoamGrid<1, 3> >::read(dune_foamgrid_path + "line1d3d.msh", /*verbose*/ true, false ) );
    {
        Dune::VTKWriter<typename FoamGrid<1, 3>::LeafGridView > writer(grid1d->leafGridView(), VTK::nonconforming);
        writer.write("1d_refined0");
    }
    gridcheck(*grid1d);
    Dune::gridinfo(*grid1d);

    // check grid adaptation interface
    checkAdaptation(*grid1d);

    // refine single element
    grid1d->mark(1, *grid1d->leafGridView().begin<0>());
    grid1d->preAdapt();
    grid1d->adapt();
    grid1d->postAdapt();

    // check grid after adaptive refinement
    checkAdaptRefinement(*grid1d);

    // write vtk
    {
        Dune::VTKWriter<typename FoamGrid<1, 3>::LeafGridView > writer(grid1d->leafGridView(), VTK::nonconforming);
        writer.write("1d_refined-l");
    }

    Dune::gridinfo(*grid1d);
    gridcheck(*grid1d);
    checkGeometryInFather(*grid1d);
    checkIntersectionIterator(*grid1d);
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
 catch (Exception e) {

    std::cout << e << std::endl;
    return 1;
 }
