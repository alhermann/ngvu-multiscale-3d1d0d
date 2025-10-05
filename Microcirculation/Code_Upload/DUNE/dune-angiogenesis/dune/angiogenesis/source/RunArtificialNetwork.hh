#include "AddFineScaleVessels.hh"
#include "AddLargeVessels.hh"
#include "ConnectVessels.hh"
#include "Plots.hh"
#include "RemoveDeadEnds.hh"
#include "SolveFluidProblem.hh"
#include "createclouds.hh"

class RunArtificialNetwork {
   public:
    typedef double RF;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<RF, 1, 1>> Matrix;
    typedef Dune::BlockVector<Dune::FieldVector<RF, 1>> Vector;
    typedef Dune::MultiTypeBlockVector<Vector, Vector> BlockVector;
    typedef Dune::MultiTypeBlockVector<Matrix, Matrix> Matrix_Row;
    typedef Dune::MultiTypeBlockMatrix<Matrix_Row, Matrix_Row> BlockMatrix;
    typedef std::vector<std::set<int>> SimplePattern;
    typedef Dune::MultiTypeBlockVector<SimplePattern, SimplePattern> Pattern_Row;
    typedef Dune::MultiTypeBlockMatrix<Pattern_Row, Pattern_Row> BlockPattern;
    typedef Dune::FoamGrid<1, 3> Grid;
    typedef Grid::LeafGridView GV;
    typedef Dune::DynamicMatrix<RF> MatrixDType;
    typedef Dune::DynamicMatrix<int> MatrixIType;
    typedef Dune::DynamicVector<int> VectorIType;

    RunArtificialNetwork() : refinement_level_3d(0), refinement_level_1d(0) {}
    virtual ~RunArtificialNetwork() {}

    void readParameters() {
        Dune::ParameterTree config_parser;
        const std::string config_filename("parameters.ini");
        Dune::ParameterTreeParser::readINITree(config_filename, config_parser);

        // 3D parameters
        refinement_level_3d = (int)config_parser.get<int>("3D.refinement_level");
        K_t = config_parser.get<RF>("3D.K_t");
        D_t = config_parser.get<RF>("3D.D_t");
        mu_int = config_parser.get<RF>("3D.mu_int");
        rho_int = config_parser.get<RF>("3D.rho_int");
        pi_int = config_parser.get<RF>("3D.pi_int");
        mO2 = config_parser.get<RF>("3D.mO2");
        cO2 = config_parser.get<RF>("3D.cO2");
        nQC = config_parser.get<RF>("3D.nQC");
        NSubdomain = config_parser.get<int>("3D.NSubdomain");
        D_tglu = config_parser.get<RF>("3D.D_tglu");
        mGlu = config_parser.get<RF>("3D.mGlu");
        mGlu_half = config_parser.get<RF>("3D.mGlu_half");
        D_tpufa = config_parser.get<RF>("3D.D_tpufa");
        mPufa = config_parser.get<RF>("3D.mPufa");
        mPufa_half = config_parser.get<RF>("3D.mPufa_half");

        // 1D parameters
        refinement_level_1d = (int)config_parser.get<int>("1D.refinement_level");
        network_file = config_parser.get<std::string>("1D.network_file");
        flow_file = config_parser.get<std::string>("1D.flow_file");
        Lp = config_parser.get<RF>("1D.Lp");
        D_v = config_parser.get<RF>("1D.D_v");
        Pv = config_parser.get<RF>("1D.Pv");
        rho_bl = config_parser.get<RF>("1D.rho_bl");
        pi_p = config_parser.get<RF>("1D.pi_p");
        sigma = config_parser.get<RF>("1D.sigma");
        gamma = config_parser.get<RF>("1D.gamma");
        threshold = config_parser.get<RF>("1D.threshold");
        mean = config_parser.get<RF>("1D.mean");
        std_dev = config_parser.get<RF>("1D.stddev");
        lambda_g = config_parser.get<RF>("1D.lambda_g");
        threshold_LV = config_parser.get<RF>("1D.threshold_LV");
        threshold_PO2 = config_parser.get<RF>("1D.threshold_PO2");
        D_vglu = config_parser.get<RF>("1D.D_vglu");
        Pv_glu = config_parser.get<RF>("1D.Pv_glu");
        D_vpufa = config_parser.get<RF>("1D.D_vpufa");
        Pv_pufa = config_parser.get<RF>("1D.Pv_pufa");
        p_arterial = config_parser.get<RF>("1D.p_arterial");

        // Transport parameter
        T_fin = config_parser.get<RF>("time.T_fin");
        dt = config_parser.get<RF>("time.dt");
        dt_explicit = config_parser.get<RF>("time.dt_explicit");

        // Solver
        solver = config_parser.get<std::string>("solver.type");
        verbosity = (int)config_parser.get<int>("solver.verbosity");
        maxIter = (int)config_parser.get<int>("solver.maxIter");
        residualReduction = config_parser.get<RF>("solver.residualReduction");

        // Plot
        print_matrix = config_parser.get<bool>("output.print_matrix");
        plot_solution = config_parser.get<bool>("output.plot_solution");
        path_plot = config_parser.get<std::string>("output.path");
        path_data = config_parser.get<std::string>("output.path_data");
        output_step = (int)config_parser.get<int>("output.output_step");
    }

    void printParameters() {
        // 3D parameters
        std::cout << "Refinement level of the 3D mesh = " << refinement_level_3d << std::endl;
        std::cout << "Permeability tissue = " << K_t << " [m^2]" << std::endl;
        std::cout << "Diffusivity tissue = " << D_t << " [m^2/s]" << std::endl;
        std::cout << "Viscosity interstitial fluid = " << mu_int << " [Pa s]" << std::endl;
        std::cout << "Density interstitial fluid = " << rho_int << unit << std::endl;
        std::cout << "Osmotic pressure interstitium = " << pi_int << " [Pa]" << std::endl;
        std::cout << "Maximal oxygen demand = " << mO2 << " [Pa/s]" << std::endl;
        std::cout << "Oxygen concentration at which the reaction rate is half of mO2 = " << cO2 << " [Pa]" << std::endl;
        std::cout << "Number of quadrature points for the coupling = " << nQC << std::endl;
        std::cout << "Diffusion parameter for glucose in tissue = " << D_tglu << " [m^2/s] " << std::endl; 
        std::cout << "Maximum consumption of glucose = " << mGlu << " [mmM/(g*s)] " << std::endl; 
        std::cout << "Half maximum consuption of glucose = " << mGlu_half << " [mmM] " << std::endl; 
        std::cout << "Diffusion parameter for PUFA in cytoplasm = " << D_tpufa << " [m^2/s] " << std::endl;
        std::cout << "Maximum consumption for PUFA = " << mPufa << " [mmM/(g*s)] " << std::endl;
        std::cout << "Half maximum consumption for PUFA = " << mPufa_half << " [mmM] " << std::endl;
        std::cout << " " << std::endl;

        // 1D parameters
        std::cout << "Refinement level of the 1D mesh " << network_file << " = " << refinement_level_1d << std::endl;
        std::cout << "Influx flow from the macrocirculation file " << network_file << std::endl;
        std::cout << "Lp = " << Lp << " [m/(Pa*s)]" << std::endl;
        std::cout << "Diffusivity vessel = " << D_v << " [m^2/s]" << std::endl;
        std::cout << "Permeability of the vessel with respect to the solute = " << Pv << " [m/s]" << std::endl;
        std::cout << "Density blood = " << rho_bl << unit << std::endl;
        std::cout << "Osmotic pressure plasma = " << pi_p << " [Pa]" << std::endl;
        std::cout << "Reflection parameter = " << sigma << " [.]" << std::endl;
        std::cout << "Murray's exponent = " << gamma << " [.]" << std::endl;
        std::cout << "Diffusion parameter for glucose in vessel (!Adjustable!) = " << D_vglu << " [m^2/s] " << std::endl;
        std::cout << "Permeability of the vessel with respect to glucose = " << Pv_glu << " [m^2/s] " << std::endl;
        std::cout << "Diffusion parameter fatty acids = " << D_vpufa << " [m^2/s] " << std::endl;
        std::cout << "Permeability of the vessel with respect to PUFA = " << Pv_pufa << " [m^2/s] " << std::endl;
        std::cout << "Minimum venous capillary pressure = " << p_arterial << " [Pa] " << std::endl;
        std::cout << " " << std::endl;

        // Transport parameter
        std::cout << "Final time = " << T_fin << " [s]" << std::endl;
        std::cout << "Time step = " << dt << " [s]" << std::endl;
        std::cout << "Time step for explicit Euler = " << dt_explicit << " [s]" << std::endl;
        std::cout << " " << std::endl;

        // Plot parameters
        std::cout << "Print matrix: " << print_matrix << std::endl;
        std::cout << "Plot solution: " << plot_solution << std::endl;
        std::cout << "Plotting solution every: " << output_step << std::endl;
        std::cout << " " << std::endl;

        // Solver's parameters
        std::cout << "Solver = " << solver << std::endl;
        std::cout << "Verbosity = " << verbosity << std::endl;
        std::cout << "Maximal iteration = " << maxIter << std::endl;
        std::cout << "Residual reduction = " << residualReduction << std::endl;
    }

    void run() {
        // Import static constants '_0', etc...
        using namespace Dune::Indices;

        const int dimworld = 3;
        int numberOfVessels_old = 0;
        int numberOfVessels_new = 1;

        Dune::Timer watch, watchTotal;
        watch.reset();
        watchTotal.reset();
        RF time = 0.0;

        std::cout << "===============================================" << std::endl;
        std::cout << "Begin simulation" << std::endl;
        std::cout << "===============================================" << std::endl;
        std::cout << " " << std::endl;
        std::cout << " " << std::endl;

        std::cout << "===============================================" << std::endl;
        std::cout << "Reading parameters..." << std::endl;
        std::cout << "===============================================" << std::endl;
        std::cout << " " << std::endl;

        readParameters();

        if (rho_bl == 1.0) {
            unit.clear();
            unit = " [m^3/s]";
        }

        std::cout << " " << std::endl;
        std::cout << " " << std::endl;

        std::cout << "===============================================" << std::endl;
        std::cout << "Printing parameters of the problem..." << std::endl;
        std::cout << "===============================================" << std::endl;
        std::cout << " " << std::endl;

        printParameters();

        std::cout << " " << std::endl;
        std::cout << " " << std::endl;

        std::cout << "===============================================" << std::endl;
        std::cout << "Creating the full 1D network..." << std::endl;
        std::cout << "===============================================" << std::endl;
        std::cout << " " << std::endl;

        watch.reset();
        std::vector<GV> networks;
        std::vector<GV> prev_networks;
        std::vector<Subdomain> subdomains;

        for (int i = 0; i < NSubdomain * NSubdomain * NSubdomain; i++) {
            subdomains.push_back(Subdomain());
        }

        networks.push_back(
            generateNetwork(vertices, indices, R, G, mu, u[_0], rho_bl, refinement_level_1d, network_file, subdomains));
        int n_segments = networks.back().size(0);
        std::pair<Dune::FieldVector<RF, dimworld>, Dune::FieldVector<RF, dimworld>> pair =
            determineExtrema(networks[0], true);
        std::pair<Dune::FieldVector<RF, dimworld>, Dune::FieldVector<RF, dimworld>> pair_o =
            determineExtrema(networks[0], false);

        // Init values
        u_init[_0] = u[_0];    
        // Init initial radii
        R0 = R;

        Dune::FieldVector<RF, dimworld> L = pair.second, O = pair.first, L_o = pair_o.second, O_o = pair_o.first;
        std::cout << "Original dimensions: L = " << L_o << ", O = " << O_o
                  << ", number of nodes = " << networks[0].size(1) << ", number of segments = " << networks[0].size(0)
                  << "\n";
        std::cout << "Enlarged domain    : L = " << L << ", O = " << O << "\n";
        std::cout << "\n";

        int numberOfEdges = indices.size();
        int numberOfEdges_new1, numberOfEdges_new2 = 0;
        int number_subdomain1, number_subdomain2, index1, index2 = 0;

        for (int i = 0; i < numberOfEdges; i++) {
            index1 = indices[i][0];
            index2 = indices[i][1];

            Dune::FieldVector<RF, 3> p1 = vertices[index1];
            Dune::FieldVector<RF, 3> p2 = vertices[index2];

            number_subdomain1 = findSubomain(p1, O, L, NSubdomain);
            number_subdomain2 = findSubomain(p2, O, L, NSubdomain);

            if (number_subdomain1 == number_subdomain2) {
                subdomains[number_subdomain1].addIndex();
                subdomains[number_subdomain1].addVertex(p1);
                subdomains[number_subdomain1].addVertex(p2);
                subdomains[number_subdomain1].addRadius(R[i]);
            } else {
                numberOfEdges_new2 = subdomains[number_subdomain2].getIndices().size();

                subdomains[number_subdomain1].addIndex();
                subdomains[number_subdomain1].addVertex(p1);
                subdomains[number_subdomain1].addVertex(p2);
                subdomains[number_subdomain1].addRadius(R[i]);

                subdomains[number_subdomain2].addIndex();
                subdomains[number_subdomain2].addVertex(p1);
                subdomains[number_subdomain2].addVertex(p2);
                subdomains[number_subdomain2].addRadius(R[i]);
            }
        }

        int step = 0;

        if (plot_solution) {
            plotNetwork(networks[0], R, G, mu, u[_0], path_plot, "network", step, refinement_level_1d);
        }
        time = watch.elapsed();
        std::cout << "Creating the 1D network took " << time << " seconds" << std::endl;
        std::cout << " " << std::endl;
        std::cout << " " << std::endl;

        std::cout << "===============================================" << std::endl;
        std::cout << "Building the 3D grid that includes the network..." << std::endl;
        std::cout << "===============================================" << std::endl;
        std::cout << " " << std::endl;

        // Create an hexahedral mesh
        std::bitset<dimworld> periodic(false);
        int overlap = 0;
        Dune::array<int, dimworld> s(Dune::fill_array<int, dimworld>(1));
        typedef typename Dune::YaspGrid<dimworld, Dune::EquidistantOffsetCoordinates<RF, dimworld>> HexahedralGrid;
        typedef typename HexahedralGrid::LeafGridView HexahedralGV;

        // Original tissue portion
        HexahedralGrid originalElement(O_o, L_o, s, periodic, overlap);
        const HexahedralGV &originalElementGV = originalElement.leafGridView();

        // Extended tissue portion
        HexahedralGrid finiteVolumeGrid(O, L, s, periodic, overlap);
        finiteVolumeGrid.globalRefine(refinement_level_3d);
        const HexahedralGV &finiteVolumeGV = finiteVolumeGrid.leafGridView();

        std::cout << "Number of elements in the 3D finite volume mesh = " << finiteVolumeGV.size(0) << "\n";
        std::cout << "Number of vertices in the 3D finite volume mesh = " << finiteVolumeGV.size(dimworld) << "\n";
        std::cout << " " << std::endl;
        std::cout << " " << std::endl;

        std::cout << "===============================================" << std::endl;
        std::cout << "Building the bounding box tree..." << std::endl;
        std::cout << "===============================================" << std::endl;
        std::cout << " " << std::endl;

        watch.reset();
        BoundingBoxTree<HexahedralGV> finiteVolumeTree(finiteVolumeGV);
        BoundingBoxTree<HexahedralGV> tree_original(originalElementGV);
        time = watch.elapsed();
        std::cout << "Determine the tree took " << time << " seconds" << std::endl;
        std::cout << " " << std::endl;
        std::cout << " " << std::endl;

        std::cout << "===============================================" << std::endl;
        std::cout << "Reading macrocirculation flow..." << std::endl;
        std::cout << "===============================================" << std::endl;
        std::cout << " " << std::endl;

        // Load the flow data into pica_flow
        loadFlowData(flow_file, pica_flow);

        // Access the flow data (example)
        for (size_t i = 0; i < pica_flow.size(); ++i) {
            std::cout << "Flow at index " << i << ": " << pica_flow[i][0] << std::endl;
        }

        std::cout << "===============================================" << std::endl;
        std::cout << "Reading neurons..." << std::endl;
        std::cout << "===============================================" << std::endl;
        std::cout << " " << std::endl;

        // Open the file with neuronal point positions
        std::ifstream neuronal_positions_file(path_plot + "neuronal_positions.csv");

        if (!neuronal_positions_file)
        {
            std::cerr << "Could not open the file with neuronal positions!" << std::endl;
        }

        // Check if opening was successful
        if (neuronal_positions_file.is_open())
        {
            std::string line;
            while (std::getline(neuronal_positions_file, line))
            {
                std::istringstream iss(line);
                Dune::FieldVector<RF, 3> neuron_point;

                for (int i = 0; i < dimworld; i++)
                {
                    std::string pos;
                    if (getline(iss, pos, ','))  // ',' is the delimiter
                    {
                        neuron_point[i] = std::stod(pos);  // Convert string to double
                    }
                }

                neuronal_vertices.push_back(neuron_point);
            }
            neuronal_positions_file.close();
        }

        // Store the vertices in a dynamic DUNE matrix
        neuronal_vertices_mat.resize(neuronal_vertices.size(), dimworld);
        for (size_t i = 0; i < neuronal_vertices.size(); i++)
        {
            for (int j = 0; j < dimworld; j++)
            {
                neuronal_vertices_mat[i][j] = neuronal_vertices[i][j];
            }
        }

        std::cout << "Neuronal points read. Number of neurons: " << neuronal_vertices.size() << "\n" << std::endl;

        // Write the matrix for center points
        center_points.resize(finiteVolumeGV.size(0), dimworld);

        size_t idx = 0;

        typedef typename HexahedralGV::Codim<0>::Iterator ElementIterator;
        for (ElementIterator eIt = finiteVolumeGV.template begin<0>(); 
                                eIt != finiteVolumeGV.template end<0>(); ++eIt, ++idx) 
        {
            const auto& x_i = eIt->geometry().center();

            for (int j = 0; j < dimworld; ++j) 
            {
                center_points[idx][j] = x_i[j];
            }
        }

        // Save the 3D center points to a file in the mounted coupling directory
        writeMatrixToFile(center_points, path_plot + "center_points_tissue.csv");
 
        // Set some simulation flags
        bool isStationary = false;
        bool continueToCalculate = true;

        std::cout << "===============================================" << std::endl;
        std::cout << "STARTING THE SIMULATION..." << std::endl;
        std::cout << "===============================================" << std::endl;
        std::cout << " " << std::endl;

        int iter = 0;
        size_t lvo = 0;
        size_t recal = 0;

        size_t lvo_exec = 0;
        size_t recal_exec = 0;
        t = 0;

        while (t < T_fin - 1.E-14) {
            t += dt_explicit;

            std::cout << "===============================================" << std::endl;
            std::cout << "CALCULATING SOLUTION AT STEP " << iter << std::endl;
            std::cout << "===============================================" << std::endl;
            std::cout << " " << std::endl;

            u[_0] *= 0.0;
            u[_1] *= 0.0;

            std::vector<ControlVolume> cv =
                solveFluidProblem(networks.back(), finiteVolumeGV, finiteVolumeTree, R, R0, G, u[_0], u_init[_0], pica_flow, path_plot, path_data,
                                solver, pattern, A, b, u, t, 20e-6, Lp, sigma, pi_p, pi_int, mu_int, rho_int, K_t,
                                residualReduction, nQC, dt_explicit, p_arterial, verbosity, maxIter, iter, output_step, plot_solution, print_matrix);

            std::cout << "Setting initial conditions...\n";
            if (iter == 0) {
                initialConditions(networks.back(), pO2, Glu, Pufa, u, R);
            }
            
            initializeMatrix(A_transport, pattern);

            // Initialize LHS matrix for glucose
            initializeMatrix(A_transport_glu, pattern);

            // Initialize LHS matrix for PUFA
            initializeMatrix(A_transport_pufa, pattern);

            BlockVector pO2_old = pO2, pO2_old_supp = pO2;
            
            // Blockvectors for glucose
            BlockVector Glu_old = Glu, Glu_old_supp = Glu;

            // Blockvectors for PUFA 
            BlockVector Pufa_old = Pufa, Pufa_old_supp = Pufa;                    

            std::cout << "===============================================" << std::endl;
            std::cout << "Perform transport simulation..." << std::endl;
            std::cout << "===============================================" << std::endl;
            std::cout << " " << std::endl;

            int inner_iter = 0;
            while (inner_iter < maxIter) {
                // Solve for glucose 

                std::cout << "Computing solution at step " << inner_iter << " seconds. Assembling transport matrix for pO2...\n";
                assembleMatrixTransport(networks.back(), finiteVolumeGV, cv, A_transport, u, pO2_old, G, K_t, mu_int,
                                        D_t, D_v, Pv, Lp, (1.0 - sigma) * (pi_p - pi_int), mO2, cO2, dt, t, 20e-6);
                std::cout << "Matrix assembled\n";
                updateRHS(networks.back(), finiteVolumeGV, pO2_old, u);
                std::cout << "Solving the linear system...\n";
                solve(A_transport, pO2_old, pO2, solver, path_plot, verbosity, maxIter, residualReduction,
                    print_matrix);

                std::cout << "Computing solution at step " << inner_iter << " seconds. Assembling transport matrix for glucose...\n";
                assembleMatrixTransport(networks.back(), finiteVolumeGV, cv, A_transport_glu, u, Glu_old, G, K_t, mu_int,
                                        D_tglu, D_vglu, Pv_glu, Lp, (1.0 - sigma) * (pi_p - pi_int), mGlu, mGlu_half, dt, t, 20e-6, true);
                std::cout << "Matrix assembled\n";
                updateRHS(networks.back(), finiteVolumeGV, Glu_old, u);
                std::cout << "Solving the linear system...\n";
                solve(A_transport_glu, Glu_old, Glu, solver, path_plot, verbosity, maxIter, residualReduction,
                    print_matrix);

                // Solve for PUFA 
                std::cout << "Computing solution at step " << inner_iter << " seconds. Assembling transport matrix for PUFA...\n";
                assembleMatrixTransport(networks.back(), finiteVolumeGV, cv, A_transport_pufa, u, Pufa_old, G, K_t, mu_int,
                                        D_tpufa, D_vpufa, Pv_pufa, Lp, (1.0 - sigma) * (pi_p - pi_int), mPufa, mPufa_half, dt, t, 20e-6);
                std::cout << "Matrix assembled\n";
                updateRHS(networks.back(), finiteVolumeGV, Pufa_old, u);
                std::cout << "Solving the linear system...\n";
                solve(A_transport_pufa, Pufa_old, Pufa, solver, path_plot, verbosity, maxIter, residualReduction,
                    print_matrix);

                // Test if equilibrium condition has been reached
                if (equilibriumReached(pO2[_1], pO2_old_supp[_1], Glu[_1], Glu_old_supp[_1], Pufa[_1], Pufa_old_supp[_1])) {
                    break;
                }

                pO2_old = pO2;
                pO2_old_supp = pO2;

                Glu_old = Glu;
                Glu_old_supp = Glu;

                Pufa_old = Pufa;
                Pufa_old_supp = Pufa;

                inner_iter += 1;
                resetMatrix(A_transport);
                resetMatrix(A_transport_glu);
                resetMatrix(A_transport_pufa);
                std::cout << " " << std::endl;
            }

            if ((plot_solution && iter % output_step == 0)) {
                std::cout << "Plotting solutions...\n";
                plot1DSolution(networks.back(), pO2[_0], R, path_plot, "pO2_explicit", step, iter);
                plot3DSolution(finiteVolumeGV, pO2[_1], path_plot, "pO2_explicit", step, iter);

                plot1DSolution(networks.back(), Glu[_0], R, path_plot, "Glu_explicit", step, iter, true);
                plot3DSolution(finiteVolumeGV, Glu[_1], path_plot, "Glu_explicit", step, iter, true);

                plot1DSolution(networks.back(), Pufa[_0], R, path_plot, "PUFA_explicit", step, iter, true);
                plot3DSolution(finiteVolumeGV, Pufa[_1], path_plot, "PUFA_explicit", step, iter, true);

                writeSolutionToFile(pO2[_0],  path_plot, "pO2_1D", static_cast<int>(std::round(t * 1000)));
                writeSolutionToFile(pO2[_1],  path_plot, "pO2_3D", static_cast<int>(std::round(t * 1000)));

                writeSolutionToFile(Glu[_0],  path_plot, "Glu_1D", static_cast<int>(std::round(t * 1000)));
                writeSolutionToFile(Glu[_1],  path_plot, "Glu_3D", static_cast<int>(std::round(t * 1000)));

                writeSolutionToFile(Pufa[_0], path_plot, "PUFA_1D", static_cast<int>(std::round(t * 1000)));
                writeSolutionToFile(Pufa[_1], path_plot, "PUFA_3D", static_cast<int>(std::round(t * 1000)));

                plot1DSolutionToCSV(u[_0], R, path_plot + "pressure_radius_1D", static_cast<int>(std::round(t * 1000)));
            }

            std::cout << " " << std::endl; 

            iter++;

            // Wait for the radii outputs to arrive 
            
            // Directory containing NVU_Vessel CSV files
            std::string directory = "NVU_Vessels";

            // Call the function
            waitForFileAndInterpolate(directory, t, radii_vec, R.size());

            // Debugging: Print the results
            for (size_t i = 0; i < radii_vec.size(); ++i) {
                std::cout << "Vessel " << i + 1 << ", Interpolated Radius: " << radii_vec[i][0] << std::endl;
                R[i] = 1e-6*radii_vec[i][0];
            }

            // R = radii_vec;                                                      // Assign the new radii
            recomputeGv(networks.back(), R, G, rho_bl);                         // Recompute conductances            
        }

        time = watchTotal.elapsed();
        std::cout << "Total time = " << std::fixed << time << " seconds" << std::endl;
    }

   private:
    BlockMatrix A, A_transport, A_transport_glu, A_transport_pufa;

    BlockVector b, u, u_init, pO2, Glu, Pufa, reducedVector1DPressure;

    BlockPattern pattern;

    MatrixDType center_points;              // Matrix for 3D tissue center points
    MatrixDType neuronal_vertices_mat;      // Neuronal positions matrix
    MatrixDType radii_vec;                  // Matrix for the updated radii outputs
    MatrixIType clouds;                     // Cloud matrix for neurons
    VectorIType cloud_lengths;              // Cloud lengths vector
    int max_cloud_length;                   // Variable for maximum cloud lengths

    Vector R, G, mu, R0;  // Data network as vectors: radii, conductances, viscosities

    std::vector<Dune::FieldVector<RF, 3>> vertices;
    std::vector<std::vector<unsigned int>> indices;

    std::vector<Dune::FieldVector<RF, 3>> neuronal_vertices;
    std::vector<Dune::FieldVector<RF, 1>> pica_flow;

    RF K_t, D_t, D_v, Lp, Pv, mu_int, rho_bl, rho_int, pi_int, sigma, pi_p, mO2, cO2, gamma, residualReduction,
        threshold, threshold_LV, threshold_PO2, nQC, 
        D_tglu, mGlu, mGlu_half, D_vglu, Pv_glu, D_tpufa, mPufa, mPufa_half, D_vpufa, Pv_pufa, p_arterial;
    RF flux_tissue_capillaries = 0.0, flux_capillaries_tissue = 0.0, mean = 0.0, std_dev = 0.0, lambda_g = 0.0, t = 0.0,
       dt, T_fin, dt_explicit;

    int refinement_level_3d, refinement_level_1d, verbosity, maxIter, readUntilDayN, NSubdomain, output_step;

    std::string path_plot, path_data, network_file, flow_file, solver;
    std::string unit = " [kg/s]";

    bool print_matrix, plot_solution;
};
