template <typename GV1D, typename Vector, typename RF, typename Subdomain>
auto connectNetwork(const GV1D& gv1d, Vector& R, Vector& G, Vector& mu, Vector& u, Dune::FieldVector<RF, 3> O,
                    Dune::FieldVector<RF, 3> L, RF rho_bl, std::vector<Subdomain>& subdomains, int NSubdomain) {
    typedef typename GV1D::Grid NetworkGrid;
    const int dimworld = GV1D::dimensionworld;
    const int dimnetwork = GV1D::dimension;

    // The network gv1d is a reference to a constant object. We first need to create a copy of the network
    Dune::GridFactory<NetworkGrid> factory;
    Dune::GeometryType type(1);
    for (const auto& vertex : vertices(gv1d)) factory.insertVertex(vertex.geometry().center());

    for (const auto& element : elements(gv1d)) {
        std::vector<unsigned int> cornerIDs(2);
        cornerIDs[0] = gv1d.indexSet().subIndex(element, 0, dimnetwork);
        cornerIDs[1] = gv1d.indexSet().subIndex(element, 1, dimnetwork);
        factory.insertElement(type, cornerIDs);
    }
    auto grid_copy = factory.createGrid();

    // Initialize vectors for the parameters of the new elements
    std::vector<RF> R_supp, G_supp, mu_supp;
    for (int i = 0; i < R.N(); i++) {
        R_supp.push_back(R[i]);
        G_supp.push_back(G[i]);
        mu_supp.push_back(mu[i]);
    }

    std::random_device rd;
    std::mt19937 generator(rd());

    // We try first to connect two points of the network, that are not too far away from each other
    for (const auto& element : elements(grid_copy->leafGridView())) {
        for (const auto& intersection : intersections(grid_copy->leafGridView(), element)) {
            if (intersection.boundary()) {
                Dune::FieldVector<RF, dimworld> vertex_coord = intersection.geometry().center();
                unsigned int index =
                    grid_copy->leafGridView().indexSet().subIndex(element, intersection.indexInInside(), dimnetwork);
                RF radius = R[grid_copy->leafGridView().indexSet().index(element)];

                Dune::FieldVector<RF, dimworld> vb = element.geometry().corner(1 - intersection.indexInInside());
                vb -= vertex_coord;

                for (const auto& element_n : elements(grid_copy->leafGridView())) {
                    for (const auto& intersection_n : intersections(grid_copy->leafGridView(), element_n)) {
                        unsigned int index_n = grid_copy->leafGridView().indexSet().subIndex(
                            element_n, intersection_n.indexInInside(), dimnetwork);
                        if (index == index_n || index < index_n) {
                            continue;
                        }

                        Dune::FieldVector<RF, dimworld> vertex_coord_n = intersection_n.geometry().center();
                        Dune::FieldVector<RF, dimworld> diff_vec = vertex_coord - vertex_coord_n;

                        int subDomain1 = findSubomain(vertex_coord, O, L, NSubdomain);
                        int subDomain2 = findSubomain(vertex_coord_n, O, L, NSubdomain);

                        RF dist = diff_vec.two_norm(), pressureDiff = (u[index] - u[index_n]) / u[index],
                           radius_n = R[grid_copy->leafGridView().indexSet().index(element_n)],
                           angle = std::acos(dot(vb, diff_vec) / (vb.two_norm() * diff_vec.two_norm()));

                        RF upperBound = 1.0e-5;

                        // std::normal_distribution<> normal_upper_bound_distribution(6.0e-5, 1.0e-5);
                        std::normal_distribution<> normal_upper_bound_distribution(1.0e-5, 2.0e-6);
                        upperBound = normal_upper_bound_distribution(generator);

                        if (dist < upperBound && std::abs(angle) < M_PI / 3.0 && radius < 4.5e-6 && radius_n < 4.5e-6 &&
                            pressureDiff > 0.004) {
                            bool is_colliding = false;

                            if (subDomain1 == subDomain2) {
                                is_colliding =
                                    testCollision(subdomains, subDomain1, vertex_coord_n, diff_vec, radius_n);
                            } else {
                                is_colliding =
                                    testCollision(subdomains, subDomain1, vertex_coord_n, diff_vec, radius_n);

                                if (!is_colliding) {
                                    is_colliding =
                                        testCollision(subdomains, subDomain2, vertex_coord_n, diff_vec, radius_n);

                                    int eIdx1D = grid_copy->leafGridView().indexSet().index(element);
                                    int eIdx1D_n = grid_copy->leafGridView().indexSet().index(element_n);
                                    RF new_R = (R[eIdx1D] + R[eIdx1D_n]) / 2.0;

                                    subdomains[subDomain2].addRadius(new_R);
                                    subdomains[subDomain2].addIndex();

                                    subdomains[subDomain2].addVertex(vertex_coord_n);
                                    subdomains[subDomain2].addVertex(vertex_coord);
                                }
                            }

                            if (!is_colliding) {
                                // Insert the new element
                                (*grid_copy).insertElement(Dune::GeometryType(1), {index, index_n});

                                // Insert new parameter values for the new segment
                                int eIdx1D = grid_copy->leafGridView().indexSet().index(element);
                                int eIdx1D_n = grid_copy->leafGridView().indexSet().index(element_n);
                                RF new_R = (R[eIdx1D] + R[eIdx1D_n]) / 2.0, new_mu = computeInVivoViscosity(new_R);
                                R_supp.push_back(new_R);
                                RF G_new = rho_bl * M_PI * std::pow(new_R, 4.0) / (8.0 * new_mu * dist);
                                G_supp.push_back(G_new);
                                mu_supp.push_back(new_mu);

                                subdomains[subDomain1].addRadius(new_R);
                                subdomains[subDomain1].addIndex();

                                subdomains[subDomain1].addVertex(vertex_coord_n);
                                subdomains[subDomain1].addVertex(vertex_coord);

                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    grid_copy->preGrow();

    // check mightVanish
    for (const auto& element : elements(grid_copy->levelGridView(0))) {
        checkHierarchy(element);
    }

    // Actually grow the network
    bool newElementGenerated = grid_copy->grow();
    if (!newElementGenerated) {
        return gv1d;
    }

    // Close growth process
    grid_copy->postGrow();
    auto gridView_copy = grid_copy->leafGridView();

    // Update parameter vectors
    G.resize(gridView_copy.size(0), 0.0);
    R.resize(gridView_copy.size(0), 0.0);
    mu.resize(gridView_copy.size(0), 0.0);
    for (int i = 0; i < R.N(); i++) {
        R[i] = R_supp[i];
        G[i] = G_supp[i];
        mu[i] = mu_supp[i];
    }

    std::cout << "At least two points in the network have been connected\n";
    std::cout << "The original network had " << gv1d.size(0) << " segments and " << gv1d.size(1) << " nodes\n";
    std::cout << "The new network has now  " << gridView_copy.size(0) << " segments and " << gridView_copy.size(1)
              << " nodes\n";

    return gridView_copy;
}
