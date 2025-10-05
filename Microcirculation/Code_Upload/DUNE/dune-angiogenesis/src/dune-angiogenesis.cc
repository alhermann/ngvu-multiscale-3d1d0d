#include <dune/angiogenesis/angiogenesis.hh>

int main(int argc, char **argv) {
    try {
        // Maybe initialize MPI
        Dune::MPIHelper &helper = Dune::MPIHelper::instance(argc, argv);

        std::cout << " " << std::endl;
        std::cout << " " << std::endl;
        std::cout << "===============================================" << std::endl;
        std::cout << "Start dune-angiogenesis." << std::endl;
        std::cout << "===============================================" << std::endl;
        std::cout << " " << std::endl;
        std::cout << " " << std::endl;

        // Get the rank of the current process 
        int myrank = helper.rank();

        // Get the size of the world communicator
        int mpiSize = helper.size();

        // Output basic MPI information
        if (mpiSize > 1)
        {
            std::cout << "MPI is active. Running a total of " << mpiSize << " processes." << std::endl;
        }

        // Simple argument parsing
        std::string scenario; // first non-flag argument is scenario
        for (int i = 1; i < argc; i++) {
            std::string arg = argv[i];
            if (scenario.empty())
                scenario = arg;
            else
                std::cout << "Unrecognized argument: " << arg << std::endl;
        }

        // If scenario is given:
        if (!scenario.empty()) {
            DuneAngiogenesis duneAngiogenesis;
            duneAngiogenesis.run(scenario);
        } else {
            std::cout << "No Scenario specified!\n";
        }

        return 0;

    } catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    } catch (...) {
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
