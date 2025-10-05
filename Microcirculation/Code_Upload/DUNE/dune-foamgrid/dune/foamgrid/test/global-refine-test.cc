// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#include <config.h>

#include "make2din3dgrid.hh"
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/grid/test/gridcheck.hh>
#include <dune/grid/test/checkintersectionit.hh>
#include <dune/grid/test/checkgeometryinfather.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

int main (int argc, char *argv[]) try
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    const std::string dune_foamgrid_path = std::string(DUNE_FOAMGRID_EXAMPLE_GRIDS_PATH) + "gmsh/";

    std::cout << "Checking FoamGrid<2, 2> (2d in 2d grid)" << std::endl;
    Dune::FieldVector<double,2> lower = {0,0};
    Dune::FieldVector<double,2> upper = {1,1};
    std::array<unsigned int,2> elements = {{1,1}};
    std::shared_ptr<FoamGrid<2, 2> > grid2d = StructuredGridFactory<FoamGrid<2,2> >::createSimplexGrid(lower,upper,elements);

    {
        Dune::VTKWriter<typename FoamGrid<2, 2>::LeafGridView >
            writer(grid2d->leafGridView(), VTK::nonconforming);
        writer.write("2d_refined0");
    }
    gridcheck(*grid2d);

    // call globalRefine with 0 (no refinement)
    // who knows who does that
    grid2d->globalRefine(0);
    gridcheck(*grid2d);

    // refine once
    grid2d->globalRefine(1);
    Dune::gridinfo(*grid2d);
    {
        Dune::VTKWriter<typename FoamGrid<2, 2>::LeafGridView >
            writer(grid2d->leafGridView(), VTK::nonconforming);
        writer.write("2d_refined1");
    }
    checkGeometryInFather(*grid2d);

    // coarsen once
    grid2d->globalRefine(-1);

    {
        Dune::VTKWriter<typename FoamGrid<2, 2>::LeafGridView >
            writer(grid2d->leafGridView(), VTK::nonconforming);
        writer.write("2d_refined-1");
    }
    checkGeometryInFather(*grid2d);
    Dune::gridinfo(*grid2d);

    // refine twice
    grid2d->globalRefine(2);
    {
        Dune::VTKWriter<typename FoamGrid<2, 2>::LeafGridView >
            writer(grid2d->leafGridView(), VTK::nonconforming);
        writer.write("2d_refined2");
    }
    Dune::gridinfo(*grid2d);
    gridcheck(*grid2d);

    // refine three times
    grid2d->globalRefine(3);
    {
        Dune::VTKWriter<typename FoamGrid<2, 2>::LeafGridView >
            writer(grid2d->leafGridView(), VTK::nonconforming);
        writer.write("2d_refined3");
    }
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

    // call globalRefine with 0 (no refinement)
    // who knows who does that
    grid1d->globalRefine(0);
    gridcheck(*grid1d);

    // refine once
    grid1d->globalRefine(); // the default is globalRefine(1)
    {
        Dune::VTKWriter<typename FoamGrid<1, 3>::LeafGridView > writer(grid1d->leafGridView(), VTK::nonconforming);
        writer.write("1d_refined1");
    }
    gridcheck(*grid1d);
    checkGeometryInFather(*grid1d);
    Dune::gridinfo(*grid1d);

    // coarsen once
    grid1d->globalRefine(-1);

    {
        Dune::VTKWriter<typename FoamGrid<1, 3>::LeafGridView > writer(grid1d->leafGridView(), VTK::nonconforming);
        writer.write("1d_refined-1");
    }
    gridcheck(*grid1d);
    checkGeometryInFather(*grid1d);
    Dune::gridinfo(*grid1d);

    // refine twice
    grid1d->globalRefine(2);
    {
        Dune::VTKWriter<typename FoamGrid<1, 3>::LeafGridView > writer(grid1d->leafGridView(), VTK::nonconforming);
        writer.write("1d_refined2");
    }
    gridcheck(*grid1d);
    Dune::gridinfo(*grid1d);

    // refine three times
    grid1d->globalRefine(3);
    {
        Dune::VTKWriter<typename FoamGrid<1, 3>::LeafGridView > writer(grid1d->leafGridView(), VTK::nonconforming);
        writer.write("1d_refined3");
    }
    gridcheck(*grid1d);
    checkIntersectionIterator(*grid1d);
    Dune::gridinfo(*grid1d);

}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
 catch (Exception e) {

    std::cout << e << std::endl;
    return 1;
 }
