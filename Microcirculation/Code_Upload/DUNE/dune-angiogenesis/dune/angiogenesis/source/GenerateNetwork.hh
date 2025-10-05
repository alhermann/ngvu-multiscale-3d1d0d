#include "SupportFunctions.hh"

template <typename Vector, typename RF, typename Subdomain>
auto generateNetwork(std::vector<Dune::FieldVector<RF, 3>>& vertices, std::vector<std::vector<unsigned int>>& indices,
                     Vector& R, Vector& G, Vector& mu, Vector& P, RF rho_bl, int refinement_level_1d,
                     const std::string network_file, std::vector<Subdomain>& subdomains) {
    const int dimnetwork = 1;
    const int dimworld = 3;
    typedef Dune::FoamGrid<dimnetwork, dimworld> Network;
    typedef Network::LeafGridView NetworkGV;
    Dune::GridFactory<Network> factory;
    Dune::GeometryType type(1);

    std::string dgf_filename(network_file);

    std::ostringstream out_dgf_filename;
    out_dgf_filename << dgf_filename.c_str();
    std::string file_dgf;
    file_dgf = out_dgf_filename.str().c_str();
    std::ifstream in_dgf_filename;
    in_dgf_filename.open(file_dgf.c_str());

    std::vector<RF> P_v, R_v, G_v, mu_v;
    std::string header;
    int fileSeparator = 0, i = 0;

    while (std::getline(in_dgf_filename, header, in_dgf_filename.widen('\n'))) {
        if (header == "# " || header == "#") {
            fileSeparator += 1;
            continue;
        } else if (i >= 3 && fileSeparator == 0) {
            Dune::FieldVector<RF, dimworld + 1> vertexInfo;

            // Split header
            std::istringstream input(header);
            std::vector<std::array<char, 50>> line;
            for (std::array<char, 50> a; input.getline(&a[0], 50, ' ');) line.push_back(a);
            int index = 0;
            for (auto& a : line) {
                vertexInfo[index] = std::atof(&a[0]);
                index += 1;
            }

            // Insert vertex
            Dune::FieldVector<RF, dimworld> vertex;
            vertex[0] = vertexInfo[0];
            vertex[1] = vertexInfo[1];
            vertex[2] = vertexInfo[2];

            factory.insertVertex(vertex);
            vertices.push_back(vertex);

            P_v.push_back(vertexInfo[3]);

        } else if (fileSeparator == 1) {
            // Avoid line: SIMPLEX in the dgf file
            fileSeparator += 1;
        } else if (fileSeparator == 2) {
            // Avoid line: parameters ... in the dgf file
            fileSeparator += 1;
        } else if (i >= 3 && fileSeparator == 3) {
            Dune::FieldVector<RF, dimworld + 1> segmentInfo;

            // Split header
            std::istringstream input(header);
            std::vector<std::array<char, 50>> line;
            for (std::array<char, 50> a; input.getline(&a[0], 50, ' ');) {
                line.push_back(a);
            }
            int index = 0;
            for (auto& a : line) {
                segmentInfo[index] = std::atof(&a[0]);
                index += 1;
            }

            // Insert segment
            std::vector<unsigned int> cornerIDs(2);
            cornerIDs[0] = (int)segmentInfo[0];
            cornerIDs[1] = (int)segmentInfo[1];
            RF radius = segmentInfo[2];

            indices.push_back(cornerIDs);
            factory.insertElement(type, cornerIDs);
            R_v.push_back(radius);
            RF viscosity = computeInVivoViscosity(radius);
            mu_v.push_back(viscosity);
            // The value still has to be divided by the length!!! Done once the network is generated.
            G_v.push_back(rho_bl * M_PI * std::pow(radius, 4.0) / (8.0 * viscosity));
        }

        i++;
    }

    int numberOfSegments = R_v.size(), numberOfPoints = P_v.size();
    P.resize(numberOfPoints);
    R.resize(numberOfSegments);
    mu.resize(numberOfSegments);
    G.resize(numberOfSegments);

    // Copy the vectors into Dune::FieldVector format
    for (int j = 0; j < numberOfPoints; j++) {
        P[j] = P_v[j];
    }
    for (int j = 0; j < numberOfSegments; j++) {
        R[j] = R_v[j];
        G[j] = G_v[j];
        mu[j] = mu_v[j];
    }

    // Keep the coarse network in order to extract correctly the pressure values
    auto network = factory.createGrid();
    network->globalRefine(refinement_level_1d);
    auto networkgv = network->leafGridView();

    // Adapt the vectors to the refinement
    Vector R_refined, G_refined, mu_refined, P_refined;
    R_refined.resize(networkgv.size(0), 0.0);
    G_refined.resize(networkgv.size(0), 0.0);
    mu_refined.resize(networkgv.size(0), 0.0);
    P_refined.resize(networkgv.size(1), 0.0);

    for (const auto& element : elements(networkgv)) {
        // Get index of the coarse element for setting radius, permeability and viscosity
        auto coarse_element = getCoarseElement(element);
        int idxElementCoarse = (network->levelGridView(0)).indexSet().index(coarse_element);
        R_refined[networkgv.indexSet().index(element)] = R[idxElementCoarse];
        G_refined[networkgv.indexSet().index(element)] = G[idxElementCoarse] / element.geometry().volume();
        mu_refined[networkgv.indexSet().index(element)] = mu[idxElementCoarse];

        // Get index for the boundary points for setting boundary conditions
        for (const auto& intersection : intersections(networkgv, element)) {
            if (intersection.boundary()) {
                int idxGlobalRefined = networkgv.indexSet().subIndex(element, intersection.indexInInside(), dimnetwork);

                // go up to the 0 level vertex
                auto level0element = intersection.inside();
                level0element = getCoarseElement(level0element);
                int idxGlobalCoarse = -1;
                for (const auto& intersection_coarse : intersections(network->levelGridView(0), level0element)) {
                    if (intersection_coarse.boundary())
                        if (Dune::FloatCmp::eq(
                                intersection_coarse.centerUnitOuterNormal() * intersection.centerUnitOuterNormal(), 1.0,
                                1.5e-7))
                            idxGlobalCoarse =
                                (network->levelGridView(0))
                                    .indexSet()
                                    .subIndex(level0element, intersection_coarse.indexInInside(), dimnetwork);
                }
                P_refined[idxGlobalRefined] = P[idxGlobalCoarse];
            }
        }
    }

    // Update vectors
    R.resize(networkgv.size(0), 0.0);
    G.resize(networkgv.size(0), 0.0);
    mu.resize(networkgv.size(0), 0.0);
    P.resize(networkgv.size(1), 0.0);

    R = R_refined;
    G = G_refined;
    mu = mu_refined;
    P = P_refined;

    return networkgv;
}

template <typename GV1D, typename BlockVector, typename Vector>
void initialConditionsFromPrevious(GV1D& previousNetwork, GV1D& currentNetwork, 
                                   BlockVector& pO2, BlockVector& Glu, BlockVector& Pufa, BlockVector& u, Vector& R) {
    using namespace Dune::Indices;
    const int dimnetwork = GV1D::dimension;
    typedef double RF;

    // Step 1: Make copies of the vectors associated with the previous network
    BlockVector pO2Prev = pO2;
    BlockVector GluPrev = Glu;
    BlockVector PufaPrev = Pufa;

    // Step 2: Resize and initialize vectors for the current network
    pO2[_0].resize(u[_0].N(), 0.0);
    pO2[_1].resize(u[_1].N(), 0.0);
    Glu[_0].resize(u[_0].N(), 0.0);
    Glu[_1].resize(u[_1].N(), 0.0);
    Pufa[_0].resize(u[_0].N(), 0.0);
    Pufa[_1].resize(u[_1].N(), 0.0);

    // Step 3: Loop over current network nodes
    for (const auto& element1d : elements(currentNetwork)) {
        auto geometryCurrent = element1d.geometry();
        for (size_t vertexIndex = 0; vertexIndex < geometryCurrent.corners(); ++vertexIndex) {
            Dune::FieldVector<RF, 3> coordCurrent = geometryCurrent.corner(vertexIndex);
            int currNodeIndex = currentNetwork.indexSet().subIndex(element1d, vertexIndex, dimnetwork);

            bool found = false; // Flag to check if the current node exists in the previous network

            // Step 4: Loop over previous network nodes to compare coordinates
            for (const auto& prevElement1d : elements(previousNetwork)) {
                auto geometryPrevious = prevElement1d.geometry();
                for (size_t prevVertexIndex = 0; prevVertexIndex < geometryPrevious.corners(); ++prevVertexIndex) {
                    Dune::FieldVector<RF, 3> coordPrevious = geometryPrevious.corner(prevVertexIndex);
                    int prevNodeIndex = previousNetwork.indexSet().subIndex(prevElement1d, prevVertexIndex, dimnetwork);

                    // Compare coordinates
                    if (coordCurrent == coordPrevious) {
                        // If the coordinates match, copy the values from the previous network's copies
                        pO2[_0][currNodeIndex] = pO2Prev[_0][prevNodeIndex];
                        Glu[_0][currNodeIndex] = GluPrev[_0][prevNodeIndex];
                        Pufa[_0][currNodeIndex] = PufaPrev[_0][prevNodeIndex];
                        found = true;
                        break; // Exit the inner loop
                    }
                }
                if (found) break; // Exit the outer loop
            }

            // Step 5: If the node is not found in the previous network, use default initialization
            if (!found) {
                RF avg = 0.0; // Compute average velocity 
                for (const auto& element1d : elements(currentNetwork))
                    avg += 0.5 * (u[_0][currentNetwork.indexSet().subIndex(element1d, 0, dimnetwork)] +
                                    u[_0][currentNetwork.indexSet().subIndex(element1d, 1, dimnetwork)]);
                avg /= currentNetwork.size(0);
                // Check intersections at this vertex
                for (const auto& intersection : intersections(currentNetwork, element1d)) {
                    if (intersection.boundary() && intersection.indexInInside() == (int)vertexIndex) {
                        // This vertex corresponds to a boundary intersection
                        int index0 = currentNetwork.indexSet().subIndex(element1d, intersection.indexInInside(), dimnetwork);
                        int index1 = currentNetwork.indexSet().subIndex(element1d, std::abs(1 - intersection.indexInInside()), dimnetwork);

                        int elementNumber = currentNetwork.indexSet().index(element1d);

                        RF velocity = u[_0][index0] - u[_0][index1];

                        /*if( velocity>=0.0 && u[_0][index0]>avg )
                        { // Arterial inflow node
                            pO2[_0][currNodeIndex ] = 75.0 * 133.32239; // [Pa], value after Pries et. al. (2013)
                            Glu[_0][currNodeIndex] =  6.0;                // Normal values for D-Glucose conc - 4 -- 6 [mM]  
                            Pufa[_0][currNodeIndex] = 0.5;               // 0.1-0.6 mM in males and 0.1-0.45 mM in females (other ref. 0.1-1 mM)
                        }
                        if( velocity>=0.0 && u[_0][index0]<avg )
                        { // Venous inflow node
                            pO2[_0][currNodeIndex ] = 38.0 * 133.32239; // [Pa], value after Pries et. al. (2013)
                            Glu[_0][currNodeIndex] = 2.5;
                            Pufa[_0][currNodeIndex] = 0.35;             // 0.35-0.5 mM venous end (males & females) 
                        }*/
                        if (R[elementNumber] > 4.4e-6)
                        {
                            // Arterial inflow node
                            pO2[_0][index0] = 75.0 * 133.32239; // [Pa], value after Pries et. al. (2013)
                            Glu[_0][index0] =  6.0;                // Normal values for D-Glucose conc - 4 -- 6 [mM]  
                            Pufa[_0][index0] = 0.5;               // 0.1-0.6 mM in males and 0.1-0.45 mM in females (other ref. 0.1-1 mM)
                        }
                        else if (R[elementNumber] < 4.0e-6 && R[elementNumber] > 2.9e-6){
                            // Arterial inflow node
                            pO2[_0][index0] = 75.0 * 133.32239; // [Pa], value after Pries et. al. (2013)
                            Glu[_0][index0] =  6.0;                // Normal values for D-Glucose conc - 4 -- 6 [mM]  
                            Pufa[_0][index0] = 0.5;               // 0.1-0.6 mM in males and 0.1-0.45 mM in females (other ref. 0.1-1 mM)
                        }
                        else if (R[elementNumber] < 4.0e-6) {
                            // Venous inflow node
                            pO2[_0][index0] = 38.0 * 133.32239; // [Pa], value after Pries et. al. (2013)
                            Glu[_0][index0] = 2.5;
                            Pufa[_0][index0] = 0.35;             // 0.35-0.5 mM venous end (males & females) 
                        }
                        else {
                            // Arterial inflow node
                            pO2[_0][index0] = 75.0 * 133.32239; // [Pa], value after Pries et. al. (2013)
                            Glu[_0][index0] =  6.0;                // Normal values for D-Glucose conc - 4 -- 6 [mM]  
                            Pufa[_0][index0] = 0.5;               // 0.1-0.6 mM in males and 0.1-0.45 mM in females (other ref. 0.1-1 mM)
                        }
                    }
                }
            }
        }
    }
}


template <typename GV1D, typename BlockVector, typename Vector>
void initialConditions(GV1D& gv1d, BlockVector& pO2, BlockVector& Glu, BlockVector& Pufa, BlockVector& u, Vector& R) {
    using namespace Dune::Indices;
    const int dimnetwork = GV1D::dimension;
    typedef double RF;

    // Initialization
    pO2[_0].resize(u[_0].N(), 0.0);
    pO2[_0] = 0.0;
    pO2[_1].resize(u[_1].N(), 0.0);
    pO2[_1] = 0.0 * 133.32239;  // [Pa]

    Glu[_0].resize(u[_0].N(), 0.0);
    Glu[_0] = 0.0;
    Glu[_1].resize(u[_1].N(), 0.0);
    Glu[_1] = 0.0;

    Pufa[_0].resize(u[_0].N(), 0.0);
    Pufa[_0] = 0.0;
    Pufa[_1].resize(u[_1].N(), 0.0);
    Pufa[_1] = 0.0;

    // double Ho2 = 0.031; //Henry's law solubility coefficient: mL of O2 / (L of blood * mmHg)
    RF avg = 0.0;

    for (const auto& element1d : elements(gv1d))
        avg += 0.5 * (u[_0][gv1d.indexSet().subIndex(element1d, 0, dimnetwork)] +
                      u[_0][gv1d.indexSet().subIndex(element1d, 1, dimnetwork)]);
    avg /= gv1d.size(0);

    // Boundary conditions on arteries and veins 
    for (const auto& element1d : elements(gv1d)) {
        for (const auto& intersection : intersections(gv1d, element1d)) {
            if (intersection.boundary()) {
                int index0 = gv1d.indexSet().subIndex(element1d, intersection.indexInInside(), dimnetwork),
                    index1 =
                        gv1d.indexSet().subIndex(element1d, std::abs(1 - intersection.indexInInside()), dimnetwork);
            
                int elementNumber = gv1d.indexSet().index(element1d);

                RF velocity = u[_0][index0] - u[_0][index1];

                if (R[elementNumber] > 4.4e-6)
                {
                    // Arterial inflow node
                    pO2[_0][index0] = 75.0 * 133.32239; // [Pa], value after Pries et. al. (2013)
                    Glu[_0][index0] =  6.0;                // Normal values for D-Glucose conc - 4 -- 6 [mM]  
                    Pufa[_0][index0] = 0.5;               // 0.1-0.6 mM in males and 0.1-0.45 mM in females (other ref. 0.1-1 mM)
                }
                else if (R[elementNumber] < 4.0e-6 && R[elementNumber] > 2.9e-6){
                    // Arterial inflow node
                    pO2[_0][index0] = 75.0 * 133.32239; // [Pa], value after Pries et. al. (2013)
                    Glu[_0][index0] =  6.0;                // Normal values for D-Glucose conc - 4 -- 6 [mM]  
                    Pufa[_0][index0] = 0.5;               // 0.1-0.6 mM in males and 0.1-0.45 mM in females (other ref. 0.1-1 mM)
                }
                else if (R[elementNumber] < 4.0e-6) {
                    // Venous inflow node
                    pO2[_0][index0] = 38.0 * 133.32239; // [Pa], value after Pries et. al. (2013)
                    Glu[_0][index0] = 2.5;
                    Pufa[_0][index0] = 0.35;             // 0.35-0.5 mM venous end (males & females) 
                }
                else {
                    // Arterial inflow node
                    pO2[_0][index0] = 75.0 * 133.32239; // [Pa], value after Pries et. al. (2013)
                    Glu[_0][index0] =  6.0;                // Normal values for D-Glucose conc - 4 -- 6 [mM]  
                    Pufa[_0][index0] = 0.5;               // 0.1-0.6 mM in males and 0.1-0.45 mM in females (other ref. 0.1-1 mM)
                }
  
                
                /*if( velocity>=0.0 && u[_0][index0]>avg )
                { // Arterial inflow node
                    pO2[_0][index0] = 75.0 * 133.32239; // [Pa], value after Pries et. al. (2013)
                    Glu[_0][index0] =  6.0;                // Normal values for D-Glucose conc - 4 -- 6 [mM]  
                    Pufa[_0][index0] = 0.5;               // 0.1-0.6 mM in males and 0.1-0.45 mM in females (other ref. 0.1-1 mM)
                }
                if( velocity>=0.0 && u[_0][index0]<avg )
                { // Venous inflow node
                    pO2[_0][index0] = 38.0 * 133.32239; // [Pa], value after Pries et. al. (2013)
                    Glu[_0][index0] = 2.5;
                    Pufa[_0][index0] = 0.35;             // 0.35-0.5 mM venous end (males & females) 
                }*/
                /*if (u[_0][index0] > avg) {                // Arterial inflow node
                    pO2[_0][index0] = 75.0 * 133.32239;   // [Pa], value after Pries et. al. (2013)
                    Glu[_0][index0] =  6.0;                // Normal values for D-Glucose conc - 4 -- 6 [mM]  
                    Pufa[_0][index0] = 0.5;               // 0.1-0.6 mM in males and 0.1-0.45 mM in females (other ref. 0.1-1 mM)
                }
                if (u[_0][index0] < avg) {               // Venous inflow node
                    pO2[_0][index0] = 38.0 * 133.32239;  // [Pa], value after Pries et. al. (2013)
                    Glu[_0][index0] = 2.5;
                    Pufa[_0][index0] = 0.35;             // 0.35-0.5 mM venous end (males & females) 
                }*/
            }
        }
    }
}

template <typename GV1D, typename BlockVector, typename Vector>
void boundaryConditions(GV1D& gv1d, BlockVector& pO2, BlockVector& Glu, BlockVector& Pufa, BlockVector& u, const Vector& R) {
    using namespace Dune::Indices;
    const int dimnetwork = GV1D::dimension;
    typedef double RF;

    // double Ho2 = 0.031; //Henry's law solubility coefficient: mL of O2 / (L of blood * mmHg)
    RF avg = 0.0;

    for (const auto& element1d : elements(gv1d))
        avg += 0.5 * (u[_0][gv1d.indexSet().subIndex(element1d, 0, dimnetwork)] +
                      u[_0][gv1d.indexSet().subIndex(element1d, 1, dimnetwork)]);
    avg /= gv1d.size(0);

    // Boundary conditions on arteries and veins 
    for (const auto& element1d : elements(gv1d)) {
        for (const auto& intersection : intersections(gv1d, element1d)) {
            if (intersection.boundary()) {
                int index0 = gv1d.indexSet().subIndex(element1d, intersection.indexInInside(), dimnetwork),
                    index1 =
                        gv1d.indexSet().subIndex(element1d, std::abs(1 - intersection.indexInInside()), dimnetwork);

                RF velocity = u[_0][index0] - u[_0][index1];
                if (u[_0][index0] > avg) {                // Arterial inflow node
                    pO2[_0][index0] = 75.0 * 133.32239;   // [Pa], value after Pries et. al. (2013)
                    Glu[_0][index0] = 6.0;
                    Pufa[_0][index0] = 0.5;               // 0.1-0.6 mM in males and 0.1-0.45 mM in females (other ref. 0.1-1 mM)
                    }
                if (u[_0][index0] < avg) {               // Venous inflow node
                    pO2[_0][index0] = 38.0 * 133.32239;  // [Pa], value after Pries et. al. (2013)
                    Glu[_0][index0] = 2.5;
                    Pufa[_0][index0] = 0.35;             // 0.35-0.5 mM venous end (males & females) 
                    }
            }
        }
    }
}