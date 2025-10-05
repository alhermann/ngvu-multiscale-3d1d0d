// Function to get the periodic flow value for a given time t
template <typename Vector, typename RF>
RF getPeriodicFlow(RF t, const Vector& flow, RF dt, RF period) {
    // Wrap time t into the range [0, period)
    double t_mod = std::fmod(t, period);
    if (t_mod < 0) t_mod += period; // Handle negative time

    // Find the closest index in the flow vector
    size_t idx = static_cast<size_t>(std::round(t_mod / dt));
    if (idx >= flow.size()) {
        idx = flow.size() - 1; // Ensure index does not exceed vector size
    }

    return flow[idx];
}

template <typename BlockMatrix, typename BlockPattern>
void initializeMatrix(BlockMatrix& A, BlockPattern& blockPattern) {
    using namespace Dune::Indices;
    using namespace Dune::Hybrid;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1> > Matrix;

    std::cout << "Initializing matrix...\n";

    int N1D = blockPattern[_0][_0].size();
    int N3D = blockPattern[_1][_1].size();
    Dune::FieldVector<int, 2> N;
    N[0] = N1D;
    N[1] = N3D;

    forEach(integralRange(Dune::Hybrid::size(A)), [&](const auto i) {
        forEach(integralRange(Dune::Hybrid::size(A[i])), [&](const auto j) {
            int nnz = 0;
            for (int k = 0; k < N[i]; k++) {
                nnz += blockPattern[i][j][k].size();
            }
            A[i][j].setSize(N[i], N[j], nnz);
            A[i][j].setBuildMode(Matrix::random);
            for (int k = 0; k < N[i]; k++) {
                A[i][j].setrowsize(k, blockPattern[i][j][k].size());
            }
            A[i][j].endrowsizes();
            for (int k = 0; k < N[i]; k++) {
                std::template set<int>::iterator setend = blockPattern[i][j][k].end();
                for (std::template set<int>::iterator setit = blockPattern[i][j][k].begin(); setit != setend; setit++) {
                    A[i][j].addindex(k, *setit);
                }
            }

            A[i][j].endindices();
            A[i][j] = 0.0;
        });
    });
}

template <typename BlockMatrix>
void resetMatrix(BlockMatrix& A) {
    // Resetting matrix to 0
    Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::Hybrid::size(A)), [&](const auto i) {
        Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::Hybrid::size(A[i])),
                              [&](const auto j) { A[i][j] = 0.0; });
    });
}

template <typename GV1D, typename GV3D, typename CV, typename BlockMatrix, typename BlockVector, typename Vector, typename stdVector,
          typename RF>
void assembleMatrix(const GV1D& gv1d, const GV3D& gv3d, CV cvs, BlockMatrix& A, BlockVector& b, 
                    Vector& G, Vector& R, Vector& R0, Vector& P, Vector& P_init, stdVector& pica_flow, RF K_t, RF mu_int, RF rho_int, RF Lp, RF osmoticGradient, 
                    RF t, RF distance_threshold, RF mean_pressure_value, RF p_arterial, const std::string path_plot,
                    bool print_matrix = false) {
    const int dimnetwork = GV1D::dimension;
    const int dimworld = GV3D::dimension;

    using namespace Dune::Indices;

    b[_0].resize(gv1d.size(1));
    b[_1].resize(gv3d.size(0));
    b = 0.0;

    RF outletArea = calculateOutletArea(gv1d, R, R0, P_init, mean_pressure_value);

    // Assemble standard elliptic operator in the network
    for (const auto& element1d : elements(gv1d)) {
        int index0 = gv1d.indexSet().subIndex(element1d, 0, dimnetwork),
            index1 = gv1d.indexSet().subIndex(element1d, 1, dimnetwork),
            elementNumber = gv1d.indexSet().index(element1d);

        A[_0][_0][index0][index0] += G[elementNumber];
        A[_0][_0][index0][index1] -= G[elementNumber];
        A[_0][_0][index1][index0] -= G[elementNumber];
        A[_0][_0][index1][index1] += G[elementNumber];

    }  // end loop over 1d-elements for standard elliptic operator

    std::cout << "Standard operator in the network assembled. Assembling tissue...\n";

    // Assemble standard elliptic operator in the tissue
    for (const auto& element3d : elements(gv3d)) {
        int indexi = gv3d.indexSet().index(element3d);
        Dune::FieldVector<RF, dimworld> x_i = element3d.geometry().center();

        for (const auto& intersection : intersections(gv3d, element3d)) {
            RF volume_face = intersection.geometry().volume();

            if (!intersection.boundary()) {
                int indexj = gv3d.indexSet().index(intersection.outside());
                Dune::FieldVector<RF, dimworld> x_j = (intersection.outside()).geometry().center();

                RF h = (x_i - x_j).two_norm();

                // Permeability K_t is homogeneous and isotropic
                RF factor = rho_int * volume_face * K_t / (h * mu_int);
                A[_1][_1][indexi][indexi] += factor;
                A[_1][_1][indexi][indexj] -= factor;

            }  // end interior face

        }  // end iterator over faces

    }  // end loop over 3d-elements for standard elliptic operator

    // ASSEMBLE COUPLING
    std::cout << "Standard operators assembled. Assembling coupling...\n";

    for (auto it = cvs.begin(); it != cvs.end(); ++it) {
        // Get the control volume
        ControlVolume cv = *it;

        int index1d = it - cvs.begin();

        std::vector<int> couplingElements = cv.get3DElementsNewCoupling();
        int numberOfElements = couplingElements.size();

        std::vector<double> weights = cv.getWeights();

        // Flux constant
        RF factor = Lp * rho_int;

        // Retrieve the center position of the 1D element (index1d)
        const int dimworld = GV1D::Grid::dimensionworld; // typically 3 for FoamGrid<1,3>
        Dune::FieldVector<RF, dimworld> x_1d;
        {
            for (const auto& element1d : elements(gv1d)) {
                int currentIndex = gv1d.indexSet().index(element1d);
                if (currentIndex == index1d) {
                    auto geometry1d = element1d.geometry();
                    x_1d = geometry1d.center(); // Center is a FieldVector<RF, dimworld>
                    break; // Exit loop once the element is found
                }
            }
        }

        for (int j = 0; j < numberOfElements; j++) {
            int index3d = couplingElements[j];

            // Assuming RF is your real type and GV3D is your 3D grid view type
            const int dimworld_3d = GV3D::Grid::dimensionworld; // Typically 3 for YaspGrid<3, ...>
            Dune::FieldVector<RF, dimworld_3d> x_3d;

            for (const auto& element3d : elements(gv3d)) {
                int currentIndex3d = gv3d.indexSet().index(element3d);
                if (currentIndex3d == index3d) {
                    auto geometry3d = element3d.geometry();
                    x_3d = geometry3d.center(); // Center of the 3D element
                    break; // Exit once found
                }
            }

            A[_1][_0][index3d][index1d] -= factor * weights[j];
            A[_0][_0][index1d][index1d] += factor * weights[j];

            // Osmotic pressures
            b[_0][index1d] -= factor * weights[j] * (-osmoticGradient);
            b[_1][index3d] += factor * weights[j] * (-osmoticGradient);

            A[_1][_1][index3d][index3d] += factor * weights[j];
            A[_0][_1][index1d][index3d] -= factor * weights[j];
        }

    }  // end loop over control volumes

    // Example: Get periodic flow values for various times
    RF flowValue = getPeriodicFlow(t, pica_flow, 0.01, 1.0);
    std::cout << "Flow at time " << t << ": " << flowValue << std::endl;

    // Get full influx
    RF totalFlux = 0.0;
    for (const auto& element1d : elements(gv1d)) {
        int index0 = gv1d.indexSet().subIndex(element1d, 0, dimnetwork),
            index1 = gv1d.indexSet().subIndex(element1d, 1, dimnetwork),
            elementNumber = gv1d.indexSet().index(element1d);

        for (const auto& intersection : intersections(gv1d, element1d)) {
            if (intersection.boundary()) {
                int index = gv1d.indexSet().subIndex(element1d, intersection.indexInInside(), dimnetwork);

                // Sum of inflow flux
                if (R0[elementNumber] > 4.4e-6)
                {
                    totalFlux += std::pow((R[elementNumber] / 8.0e-4), 3) * flowValue;
                }
                else if (R0[elementNumber] < 4.0e-6 && R0[elementNumber] > 2.9e-6)
                {   
                    totalFlux += std::pow((R[elementNumber] / 8.0e-4), 3) * flowValue; 
                }
            }
        }
    }

    std::cout << "Total Flux: " << totalFlux << std::endl;

    // Impose Dirichlet boundary conditions on the network
    for (const auto& element1d : elements(gv1d)) {
        int index0 = gv1d.indexSet().subIndex(element1d, 0, dimnetwork),
            index1 = gv1d.indexSet().subIndex(element1d, 1, dimnetwork),
            elementNumber = gv1d.indexSet().index(element1d);

        for (const auto& intersection : intersections(gv1d, element1d)) {
            if (intersection.boundary()) {
                int index = gv1d.indexSet().subIndex(element1d, intersection.indexInInside(), dimnetwork);

                if (R0[elementNumber] > 4.4e-6)
                {
                    b[_0][index] = std::pow((R[elementNumber] / 8.0e-4), 3) * flowValue; 
                }
                else if (R0[elementNumber] < 4.0e-6 && R0[elementNumber] > 2.9e-6){
                    b[_0][index] = std::pow((R[elementNumber] / 8.0e-4), 3) * flowValue;
                }
                else if (R0[elementNumber] < 4.0e-6) {
                    b[_0][index] = -totalFlux * (R[elementNumber] * R[elementNumber] * M_PI) / outletArea;                  
                }
                else {
                    b[_0][index] = std::min(static_cast<RF>(P_init[index]), p_arterial);
                        
                    A[_0][_0][index] = 0.0;
                    A[_0][_1][index] = 0.0;
                    A[_0][_0][index][index] = 1.0;
                }

            }  // end over boundary vertex

        }  // end iterator on vertices

    }  // end iterator over elements

    if (print_matrix) {
        std::cout << "Printing the matrix A and the vector b...\n";
        Dune::BCRSMatrix<Dune::FieldMatrix<RF, 1, 1> > B;
        Dune::BlockVector<Dune::FieldVector<RF, 1> > F;
        packBlockMatrix(A, B);
        packBlockVector(b, F);
        std::string filenameA = path_plot, filenameb = path_plot;

        filenameA.append("A.m");
        filenameb.append("b.m");
        Dune::writeMatrixToMatlab(B, filenameA, 18);
        Dune::writeVectorToMatlab(F, filenameb, 18);
    }
}

/*
 * CCFV method for the 3D-coupled transport equation
    \nabla \cdot \{v u - D \nabla u \} = q in \Omega
    \{v u - D \nabla u \} \cdot \nu = \beta u on \partial\Omega
    Mass-flux discretizazion for the 1D-transport equation
 */
template <typename GV1D, typename GV3D, typename CV, typename BlockMatrix, typename BlockVector, typename Vector,
          typename RF>
void assembleMatrixTransport(const GV1D& gv1d, const GV3D& gv3d, CV cvs, BlockMatrix& A, BlockVector& u,
                             BlockVector& c_old, Vector& G, RF K_t, RF mu_int, RF D_t, RF D_v, RF Pv, RF Lp,
                             RF osmoticGradient, RF mO2, RF cO2, RF dt, RF t, RF distance_threshold, bool glu_case=false) {
    const int dimnetwork = GV1D::dimension;
    const int dimworld = GV3D::dimension;

    using namespace Dune::Indices;

    // Assemble standard advection and diffusion in the network
    int N1D = gv1d.size(1);

    for (int i = 0; i < N1D; i++) {
        A[_0][_0][i][i] = 1.0;
    }

    for (const auto& element1d : elements(gv1d)) {
        int eIdx1D = gv1d.indexSet().index(element1d), index0 = gv1d.indexSet().subIndex(element1d, 0, dimnetwork),
            index1 = gv1d.indexSet().subIndex(element1d, 1, dimnetwork);

        // Advective term: calculate velocity at the node under consideration
        RF vol_cv0 = cvs[index0].getVolume(), vol_cv1 = cvs[index1].getVolume();

        RF flux01 = dt * G[eIdx1D] * (u[_0][index0] - u[_0][index1]), flux10 = -flux01;

        // Upwinding
        if (flux01 < 0.0) {
            A[_0][_0][index0][index1] += flux01 / vol_cv0;
        } else {
            A[_0][_0][index0][index0] += flux01 / vol_cv0;
        }

        if (flux10 < 0.0) {
            A[_0][_0][index1][index0] += flux10 / vol_cv1;
        } else {
            A[_0][_0][index1][index1] += flux10 / vol_cv1;
        }
        
        // Diffusion term, assumed homogeneous
        RF h = element1d.geometry().volume(), factor = dt * D_v / (h * h);
        A[_0][_0][index0][index0] += factor;
        A[_0][_0][index0][index1] -= factor;
        A[_0][_0][index1][index0] -= factor;
        A[_0][_0][index1][index1] += factor;

    }  // end loop over 1d-elements for standard advection and diffusion in the network

    RF D_t_at_x = D_t; // Default value of D_t
    RF mO2_at_x = mO2; // Default value of m02

    // Assemble standard advection and diffusion in the tissue
    for (const auto& element3d : elements(gv3d)) {
        int indexi = gv3d.indexSet().index(element3d);
        Dune::FieldVector<RF, dimworld> x_i = element3d.geometry().center();

        // Determine D_t based on proximity to the 1D network
        D_t_at_x = D_t; // Default value of D_t
        mO2_at_x = mO2; // Default value of mO2

        // Contribution mass matrix
        A[_1][_1][indexi][indexi] += std::abs(element3d.geometry().volume());

        // Oxygen metabolization (implemented semi-implicitly)
        A[_1][_1][indexi][indexi] += dt * std::abs(element3d.geometry().volume()) * mO2_at_x / (cO2 + c_old[_1][indexi]);

        for (const auto& intersection : intersections(gv3d, element3d)) {
            RF volume_face = intersection.geometry().volume();

            if (!intersection.boundary()) {
                int indexj = gv3d.indexSet().index(intersection.outside());
                Dune::FieldVector<RF, dimworld> x_j = (intersection.outside()).geometry().center();

                RF h = (x_i - x_j).two_norm();

                // Advective term: calculate velocity at face's center
                RF v = -K_t * (u[_1][indexj] - u[_1][indexi]) / (h * mu_int);
                if (v >= 0) {
                    // Upwinding
                    A[_1][_1][indexi][indexi] += dt * v * volume_face;
                } else {
                    A[_1][_1][indexi][indexj] += dt * v * volume_face;
                }

                // Diffusion term: tensor assumed homogeneous and isotropic
                RF factor = dt * D_t_at_x * volume_face / h;
                A[_1][_1][indexi][indexi] += factor;
                A[_1][_1][indexi][indexj] -= factor;

            }  // end interior face
            else {
                // Inhomogeneous Neumann conditions that depend on the concentration itself
                RF beta = 0.0;  // 1e-3;
                A[_1][_1][indexi][indexi] += dt * beta * volume_face;
            }

        }  // end iterator over faces

    }  // end loop over 3d-elements for standard advection and diffusion in the tissue
    
    // ASSEMBLE COUPLING
    for (auto it = cvs.begin(); it != cvs.end(); ++it) {
        // Get the control volume
        ControlVolume cv = *it;

        int index1d = it - cvs.begin();
        std::vector<int> couplingElements = cv.get3DElementsNewCoupling();
        int numberOfElements = couplingElements.size();

        std::vector<double> weights = cv.getWeights();
        RF vol_cv = cv.getVolume();

        // Get pressure in network at the center of the control volume
        RF p_v = u[_0][index1d];

        for (int j = 0; j < numberOfElements; j++) {
            int index3d = couplingElements[j];

            // Assuming RF is your real type and GV3D is your 3D grid view type
            const int dimworld_3d = GV3D::Grid::dimensionworld; // Typically 3 for YaspGrid<3, ...>
            Dune::FieldVector<RF, dimworld_3d> x_3d;

            for (const auto& element3d : elements(gv3d)) {
                int currentIndex3d = gv3d.indexSet().index(element3d);
                if (currentIndex3d == index3d) {
                    auto geometry3d = element3d.geometry();
                    x_3d = geometry3d.center(); // Center of the 3D element
                    break; // Exit once found
                }
            }

            // Calculate filtration
            RF filtration = Lp * (p_v - osmoticGradient);

            A[_0][_0][index1d][index1d] += dt * (0.5 * (filtration - Lp * u[_1][index3d]) + Pv) * weights[j] / vol_cv;
            A[_1][_0][index3d][index1d] -= dt * (0.5 * (filtration - Lp * u[_1][index3d]) + Pv) * weights[j];

            A[_1][_1][index3d][index3d] -= dt * (0.5 * (filtration - Lp * u[_1][index3d]) - Pv) * weights[j];
            A[_0][_1][index1d][index3d] += dt * (0.5 * (filtration - Lp * u[_1][index3d]) - Pv) * weights[j] / vol_cv;
        }
    }

    // Impose known concentrations on the network's boundaries
    for (const auto& element1d : elements(gv1d)) {
        for (const auto& intersection : intersections(gv1d, element1d)) {
            if (intersection.boundary()) {
                int index0 = gv1d.indexSet().subIndex(element1d, intersection.indexInInside(), dimnetwork),
                    index1 =
                        gv1d.indexSet().subIndex(element1d, std::abs(1 - intersection.indexInInside()), dimnetwork);

                A[_0][_0][index0] = 0.0;
                A[_0][_1][index0] = 0.0;
                A[_0][_0][index0][index0] = 1.0;

                RF velocity = u[_0][index0] - u[_0][index1];
                if (velocity < 0.0) {
                    A[_0][_0][index0][index1] = -1.0;
                }
            }
        }
    }
}

template <typename GV1D, typename GV3D, typename BlockVector>
void updateRHS(const GV1D& gv1d, const GV3D& gv3d, BlockVector& b, BlockVector& u) {
    const int dimnetwork = GV1D::dimension;
    const int dimworld = GV3D::dimension;

    using namespace Dune::Indices;

    // Multiply the right-hand side by the mass-matrix
    for (const auto& element3d : elements(gv3d)) {
        int indexi = gv3d.indexSet().index(element3d);
        b[_1][indexi] *= std::abs(element3d.geometry().volume());
    }

    for (const auto& element1d : elements(gv1d)) {
        for (const auto& intersection : intersections(gv1d, element1d)) {
            if (intersection.boundary()) {
                int index0 = gv1d.indexSet().subIndex(element1d, intersection.indexInInside(), dimnetwork),
                    index1 =
                        gv1d.indexSet().subIndex(element1d, std::abs(1 - intersection.indexInInside()), dimnetwork);

                double velocity = u[_0][index0] - u[_0][index1];

                if (velocity < 0.0) {
                    b[_0][index0] = 0.0;
                }
            }
        }

    }  // end loop over 1d-elements for standard advection and diffusion in the network
}