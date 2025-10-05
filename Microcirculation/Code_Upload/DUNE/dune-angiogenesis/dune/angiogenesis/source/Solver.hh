template <typename Matrix, typename Vector>
void solve(Matrix& A, Vector& b, Vector& u, std::string solverType, const std::string path_plot, int verbosity = 1,
           int maxIter = 1000/*500*/, double residualReduction = 1e-13, bool print_matrix = false) {
    using namespace Dune::Indices;
    using namespace Dune::Hybrid;

    int N = 0;
    forEach(integralRange(Dune::Hybrid::size(A)), [&](const auto i) { N += A[i][i].N(); });

    if (verbosity > 0)
        std::cout << "Solving linear system of dimension " << N << "x" << N << " with the " << solverType << " method"
                  << std::endl;

    typedef double RF;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<RF, 1, 1> > SparseMatrix;
    typedef Dune::BlockVector<Dune::FieldVector<RF, 1> > SparseVector;
    SparseVector U;

    Dune::Timer watchSolverTotal, watchPacking;
    watchSolverTotal.reset();

    // initialize u to some arbitrary value to avoid u being the exact solution
    forEach(integralRange(Dune::Hybrid::size(A)), [&](const auto i) {
        u[i].resize(b[i].N(), false);
        u[i] = 1.0;
    });

    if (solverType == "BlockDiagAMGBiCGSTAB") {
        auto linearSolver = std::make_shared<Dune::BlockDiagAMGBiCGSTABSolver>();
        const bool converged = linearSolver->template solve<2>(A, u, b, residualReduction, maxIter, verbosity);
        if (!converged)
            std::cout << "Solver did not converged!\n";
        else if (converged && verbosity > 0)
            std::cout << "Solver converged after " << linearSolver->result().iterations
                      << " iterations with a final reduction of " << linearSolver->result().reduction << "\n";
    } else {
        std::cout << "==================================================================" << std::endl;
        std::cout << "The solver " << solverType << " is not available!!!" << std::endl;
        std::cout << "==================================================================" << std::endl;
        std::cout << " " << std::endl;
        std::cout << " " << std::endl;
        exit(0);
    }

    if (print_matrix) {
        std::cout << "Printing the solution vector u...\n";
        packBlockVector(u, U);
        std::string filename = path_plot;
        filename.append("U.m");
        Dune::writeVectorToMatlab(U, filename, 18);
    }

    RF solverTime = watchSolverTotal.elapsed();
    std::cout << "Solving the linear system took " << solverTime << " seconds" << std::endl;
}
