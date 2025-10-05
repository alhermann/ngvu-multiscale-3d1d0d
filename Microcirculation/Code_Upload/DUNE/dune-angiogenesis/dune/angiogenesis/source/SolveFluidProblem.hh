#include <math.h>
#include <time.h>

#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "AdjacencyPattern.hh"
#include "AssembleMatrix.hh"
#include "ControlVolume.hh"
#include "Fluxes.hh"
#include "ReportData.hh"

using namespace Dune::Indices;

template <typename GV1D, typename GV3D, typename Tree, typename Vector, typename stdVector, typename BlockPattern, typename BlockMatrix,
          typename BlockVector, typename RF>
std::vector<ControlVolume> solveFluidProblem(const GV1D& gv1d, const GV3D& gv3d, Tree& tree, Vector& R, Vector& R0, Vector& G,
                                             Vector& P, Vector& P_init, stdVector &pica_flow, const std::string path_plot, const std::string path_data,
                                             std::string solver, BlockPattern& blockPattern, BlockMatrix& A,
                                             BlockVector& b, BlockVector& u, RF t, RF distance_threshold,
                                             RF Lp, RF sigma, RF pi_p, RF pi_int, RF mu_int, RF rho_int, RF K_t, RF residualReduction, RF nQC, 
                                             RF dt_explicit, RF p_arterial, int verbosity, int maxIter, int step, int output_step, bool plot_solution, bool print_matrix) {
    RF threshold = 7e-6;
    RF mean_pressure_value = 4617.2;

    Dune::Timer watch;
    watch.reset();

    std::cout << "===============================================" << std::endl;
    std::cout << "Determine coupling stencils..." << std::endl;
    std::cout << "===============================================" << std::endl;
    std::cout << " " << std::endl;

    watch.reset();
    std::vector<ControlVolume> cv = determineStencils(gv1d, gv3d, tree, /*R,*/ R0, nQC, threshold);
    RF time = watch.elapsed();

    std::cout << "Determine the stencils took " << time << " seconds" << std::endl;
    std::cout << " " << std::endl;
    std::cout << " " << std::endl;
    std::cout << "===============================================" << std::endl;
    std::cout << "Determine adjacency pattern..." << std::endl;
    std::cout << "===============================================" << std::endl;
    std::cout << " " << std::endl;

    watch.reset();
    std::cout << "Determine pattern... " << std::endl;
    determinePattern(gv1d, gv3d, cv, blockPattern);
    initializeMatrix(A, blockPattern);
    time = watch.elapsed();
    std::cout << "Determine the pattern took " << time << " seconds" << std::endl;
    std::cout << " " << std::endl;
    std::cout << " " << std::endl;

    std::cout << "===============================================" << std::endl;
    std::cout << "Assembling matrix..." << std::endl;
    std::cout << "===============================================" << std::endl;
    std::cout << " " << std::endl;

    watch.reset();
    std::cout << "Assemble matrix... " << std::endl;
    assembleMatrix(gv1d, gv3d, cv, A, b, G, R, R0, P, P_init, pica_flow, K_t, mu_int, rho_int, Lp, sigma * (pi_p - pi_int), 
                   t, distance_threshold, mean_pressure_value, p_arterial, path_plot, print_matrix);

    time = watch.elapsed();
    std::cout << "Assembling the matrix took " << time << " seconds" << std::endl;
    std::cout << " " << std::endl;
    std::cout << " " << std::endl;

    std::cout << "===============================================" << std::endl;
    std::cout << "Solving the linear system..." << std::endl;
    std::cout << "===============================================" << std::endl;
    std::cout << " " << std::endl;

    solve(A, b, u, solver, path_plot, verbosity, maxIter, residualReduction, print_matrix);

    std::cout << " " << std::endl;
    std::cout << " " << std::endl;

    writeSolutionToFile(u[_0], "u_1D", step);
    writeSolutionToFile(u[_1], "u_3D", step);

    if (plot_solution) {
        std::cout << "===============================================" << std::endl;
        std::cout << "Flow equations: VTK Output" << std::endl;
        std::cout << "===============================================" << std::endl;
        std::cout << " " << std::endl;

        plot1DSolution(gv1d, u[_0], R, path_plot, "flow", step);
        std::cout << "1D solution plotted. Saving 3D solution..." << std::endl;
        plot3DSolution(gv3d, u[_1], path_plot, "flow", step);
        std::cout << " " << std::endl;
        std::cout << " " << std::endl;
    }

    std::cout << "===============================================" << std::endl;
    std::cout << "Computing fluxes..." << std::endl;
    std::cout << "===============================================" << std::endl;
    std::cout << " " << std::endl;

    RF flux_tissue_capillaries = 0.0;
    RF flux_capillaries_tissue = 0.0;

    computeFluxes(gv1d, gv3d, cv, u[_0], u[_1], rho_int, Lp, sigma * (pi_p - pi_int), flux_capillaries_tissue,
                flux_tissue_capillaries);

    std::cout << "Total flux capillaries -> tissue =  " << std::setprecision(10) << flux_capillaries_tissue * 1e9
            << " [mu g/s]\n";
    std::cout << "Total flux tissue -> capillaries = " << std::setprecision(10) << flux_tissue_capillaries * 1e9
            << " [mu g/s]\n";

    std::cout << " " << std::endl;
    std::cout << " " << std::endl;

    reportFluxIntoTissue(path_data, flux_capillaries_tissue * 1e9);

    return cv;
}
