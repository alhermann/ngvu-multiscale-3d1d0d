template <typename GV1D, typename GV3D, typename BlockVector, typename Vector, typename RF, typename Subdomain>
auto addFineScaleVessels(const GV1D& gv1d, const GV3D gv3d, Vector& R, Vector& G, Vector& mu, BlockVector& pO2,
                         Vector& u, RF rho_bl, RF gamma, RF mean, RF std_dev, RF lambda_g, Dune::FieldVector<RF, 3> O,
                         Dune::FieldVector<RF, 3> L, std::vector<Subdomain>& subdomains, int NSubdomain,
                         std::vector<RF>& avg_local_po2, RF threshold_PO2) {
    typedef typename GV1D::Grid NetworkGrid;
    const int dimworld = GV3D::dimension;
    const int dimnetwork = GV1D::dimension;
    using namespace Dune::Indices;

    std::lognormal_distribution<> log_normal_distribution(mean, std_dev);
    std::random_device rd;
    std::mt19937 generator(rd());

    int counter_new_vessels = 0;

    // Parameters needed to find the 3D elements around a point of the network
    int level_3d = gv3d.grid().maxLevel(), N = std::pow(2.0, level_3d);
    Dune::FieldVector<RF, dimworld> origin(0.0),  // Origin point
        h(0.0);                                   // Meshsizes
    for (const auto& element : elements(gv3d)) {
        for (int i = 0; i < dimworld; ++i)
            h[i] = std::abs(element.geometry().corner(0)[i] - element.geometry().corner(7)[i]);
        origin = element.geometry().corner(0);
        break;
    }

    typedef typename GV3D::template Codim<0>::Entity Element;
    IndexToElementMap<GV3D> indexToElementMap(gv3d);
    indexToElementMap.resize(gv3d.size(0));
    for (auto&& element : elements(gv3d)) {
        unsigned int eIdx = gv3d.indexSet().index(element);
        indexToElementMap[eIdx] = element.seed();
    }

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
    std::vector<RF> R_supp, G_supp, mu_supp, u_supp, pO2_supp;
    for (int i = 0; i < R.N(); i++) {
        R_supp.push_back(R[i]);
        G_supp.push_back(G[i]);
        mu_supp.push_back(mu[i]);
    }
    for (int i = 0; i < u.N(); i++) {
        u_supp.push_back(u[i]);
        pO2_supp.push_back(pO2[_0][i]);
    }

    for (const auto& element : elements(grid_copy->leafGridView())) {
        int eIdx1D = grid_copy->leafGridView().indexSet().index(element);

        RF R_p = R[eIdx1D];

        for (const auto& intersection : intersections(grid_copy->leafGridView(), element)) {
            if (intersection.boundary()) {
                // Informations about the vertex
                Dune::FieldVector<RF, dimworld> vertex_coord = intersection.geometry().center();
                unsigned int index =
                    grid_copy->leafGridView().indexSet().subIndex(element, intersection.indexInInside(), dimnetwork);

                // Vector pointing from the interior point to the boundary point
                Dune::FieldVector<RF, dimworld> vb = element.geometry().corner(1 - intersection.indexInInside());
                vb -= vertex_coord;

                // Determine corners of the hexahedral element with center in the vertex
                Dune::FieldVector<RF, dimworld> c1 = vertex_coord;
                Dune::FieldVector<RF, dimworld> c2 = vertex_coord;
                c1 -= 3e-5;  // Move the first corner of 30 mu m (see Pries et.al. [2013])
                c2 += 3e-5;  // Move the last corner of 30 mu m (see Pries et.al. [2013])

                int i1 = (c1[0] - O[0]) / h[0], j1 = (c1[1] - O[1]) / h[1], k1 = (c1[2] - O[2]) / h[2],
                    i2 = (c2[0] - O[0]) / h[0], j2 = (c2[1] - O[1]) / h[1], k2 = (c2[2] - O[2]) / h[2];

                // Check first if the corner is outside the domain
                i1 < 0 ? i1 = 0 : i1;
                i1 >= N ? i1 = N - 1 : i1;
                i2 < 0 ? i2 = 0 : i2;
                i2 >= N ? i2 = N - 1 : i2;
                j1 < 0 ? j1 = 0 : j1;
                j1 >= N ? j1 = N - 1 : j1;
                j2 < 0 ? j2 = 0 : j2;
                j2 >= N ? j2 = N - 1 : j2;
                k1 < 0 ? k1 = 0 : k1;
                k1 >= N ? k1 = N - 1 : k1;
                k2 < 0 ? k2 = 0 : k2;
                k2 >= N ? k2 = N - 1 : k2;

                RF min_pO2 = 1e5;
                int min_n = N * N * N;
                for (int i = i1; i <= i2; ++i)
                    for (int j = j1; j <= j2; ++j)
                        for (int k = k1; k <= k2; ++k) {
                            int n = i + j * N + k * N * N;

                            Element e = indexToElementMap.entity(n);
                            Dune::FieldVector<RF, dimworld> element3D_center = e.geometry().center();
                            Dune::FieldVector<RF, dimworld> v2 = element3D_center;
                            v2 -= vertex_coord;

                            // The segment was generated in the previous step and therefore the new node was the center
                            // of the 3D-element already
                            if (v2.two_norm() < 1e-16) continue;

                            // The center of the element is inside the considered cone => we can update the values
                            RF dist = 0.0;
                            if (min_n != N * N * N)
                                dist =
                                    ((indexToElementMap.entity(min_n)).geometry().center() - vertex_coord).two_norm();
                            if (pO2[_1][n] < min_pO2 - 1e-10) {
                                min_pO2 = pO2[_1][n];
                                min_n = n;
                            } else if (std::abs(pO2[_1][n] - min_pO2) < 1e-10 && v2.two_norm() > dist) {
                                min_pO2 = pO2[_1][n];
                                min_n = n;
                            }
                        }

                // No element contained in the cone: can happen at nodes on the 3D-boundary
                if (min_n == N * N * N) continue;

                RF x = log_normal_distribution(generator);

                // Do we generate a bifurcation?
                RF prob = 0.5 + 0.5 * std::erf((std::log(x) - mean) / std::sqrt(2.0 * std_dev * std_dev));
                bool bifurcation = false;
                prob > 0.6 ? bifurcation = true : false;

                RF R_c = std::pow(2.0, -1.0 / gamma) * R_p;
                R_c > R_p ? R_c = R_p : R_c;

                std::normal_distribution<> normal_distribution(R_c, R_c / 32.0);
                RF R_1 = normal_distribution(generator);

                if (R_1 < 3.0e-6)  // 3.25e-6
                {
                    std::normal_distribution<> normal_small_R_distribution(2.75e-6, 2.5e-7);  // 3.0e-6 2.5e-7
                    R_1 = normal_small_R_distribution(generator);
                    R_1 < 2.0e-6 ? R_1 = 2.0e-6 : R_1;  // 2.5e-6
                }

                if (!bifurcation) R_1 = R_p;

                // We have now the number of the element with the lowest pO2 contained in the cone around the considered
                // vertex in the network.
                Element e_min_pO2 = indexToElementMap.entity(min_n);
                Dune::FieldVector<RF, dimworld> v1 = e_min_pO2.geometry().center(), vb_normed = vb;

                vb_normed *= lambda_g;

                // It's -vb, because vb=v_interior-v_boundary and we need the opposite vector
                v1 -= vb_normed;

                Dune::FieldVector<RF, 1> s = x * R_1;
                s < 1.0e-5 ? s = 1.0e-5 : s;

                Dune::FieldVector<RF, dimworld> v1_moved =
                    pointOnSegment(vertex_coord, e_min_pO2.geometry().center(), s);

                if (isOutside(v1_moved, O, L)) v1_moved = e_min_pO2.geometry().center();

                int subDomain1 = findSubomain(vertex_coord, O, L, NSubdomain);
                int subDomain2 = findSubomain(v1_moved, O, L, NSubdomain);

                bool is_colliding = testCollision(subdomains, subDomain1, vertex_coord, v1_moved - vertex_coord, R_1);

                if (!is_colliding && (subDomain1 != subDomain2)) {
                    is_colliding = testCollision(subdomains, subDomain2, vertex_coord, v1_moved - vertex_coord, R_1);

                    if (!is_colliding && (avg_local_po2[subDomain2] < threshold_PO2 * 133.32239)) {
                        subdomains[subDomain2].addRadius(R_1);
                        subdomains[subDomain2].addIndex();
                        subdomains[subDomain2].addVertex(vertex_coord);
                        subdomains[subDomain2].addVertex(v1_moved);
                    }
                }

                if (!is_colliding && (avg_local_po2[subDomain1] < threshold_PO2 * 133.32239)) {
                    if (intersection.geometry().center()[2] > 1.1e-5) {
                        // Add vertex to the network
                        auto newVIdx = (*grid_copy).insertVertex(v1_moved);

                        // Insert the new element
                        (*grid_copy).insertElement(Dune::GeometryType(1), {index, newVIdx});

                        // Insert new parameter values for the new segment
                        RF mu_new = computeInVivoViscosity(R_1),
                           G_new = rho_bl * M_PI * std::pow(R_1, 4.0) /
                                   (8.0 * mu_new * (vertex_coord - v1_moved).two_norm());
                        R_supp.push_back(R_1);
                        G_supp.push_back(G_new);
                        mu_supp.push_back(mu_new);
                        u_supp.push_back(u[index]);
                        pO2_supp.push_back(pO2[_0][index]);

                        subdomains[subDomain1].addRadius(R_1);
                        subdomains[subDomain1].addIndex();
                        subdomains[subDomain1].addVertex(vertex_coord);
                        subdomains[subDomain1].addVertex(v1_moved);

                        counter_new_vessels += 1;
                    }
                }

                if (bifurcation) {
                    RF R_2 = normal_distribution(generator);

                    if (R_2 < 3.0e-6) {
                        std::normal_distribution<> normal_small_R_distribution(2.75e-6, 2.5e-7);
                        R_2 = normal_small_R_distribution(generator);
                        R_2 < 2.0e-6 ? R_2 = 2.0e-6 : R_2;
                    }

                    RF phi = std::acos((std::pow(R_p, 4.0) + std::pow(R_2, 4.0) - std::pow(R_1, 4.0)) /
                                       (2.0 * std::pow(R_p, 2.0) * std::pow(R_2, 2.0)));
                    Dune::FieldVector<RF, dimworld> normal =
                        crossProduct(v1_moved - vertex_coord,
                                     element.geometry().corner(1 - intersection.indexInInside()) - vertex_coord);
                    Dune::FieldVector<RF, dimworld> v_bif = findPointOnCircleBifurcation(
                        normal, vertex_coord, (v1_moved - vertex_coord).two_norm(), phi, vb);

                    if (isOutside(v_bif, O, L)) continue;

                    if (std::isnan(v_bif[0]) || std::isnan(v_bif[1]) || std::isnan(v_bif[2])) continue;

                    int subDomain2 = findSubomain(v_bif, O, L, NSubdomain);

                    is_colliding = testCollision(subdomains, subDomain1, vertex_coord, v_bif - vertex_coord, R_2);

                    if (!is_colliding && (subDomain1 != subDomain2)) {
                        is_colliding = testCollision(subdomains, subDomain2, vertex_coord, v_bif - vertex_coord, R_2);

                        if (!is_colliding && (avg_local_po2[subDomain2] < threshold_PO2 * 133.32239)) {
                            subdomains[subDomain2].addRadius(R_2);
                            subdomains[subDomain2].addIndex();
                            subdomains[subDomain2].addVertex(vertex_coord);
                            subdomains[subDomain2].addVertex(v_bif);
                        }
                    }

                    if (!is_colliding && (avg_local_po2[subDomain1] < threshold_PO2 * 133.32239)) {
                        if (intersection.geometry().center()[2] > 1.1e-5) {
                            auto newVIdxBif = (*grid_copy).insertVertex(v_bif);

                            (*grid_copy).insertElement(Dune::GeometryType(1), {index, newVIdxBif});

                            RF mu_new_n = computeInVivoViscosity(R_2),
                               G_new_n = rho_bl * M_PI * std::pow(R_2, 4.0) /
                                         (8.0 * mu_new_n * (vertex_coord - v_bif).two_norm());
                            R_supp.push_back(R_2);
                            mu_supp.push_back(mu_new_n);
                            G_supp.push_back(G_new_n);
                            u_supp.push_back(u[index]);
                            pO2_supp.push_back(pO2[_0][index]);

                            subdomains[subDomain1].addRadius(R_2);
                            subdomains[subDomain1].addIndex();
                            subdomains[subDomain1].addVertex(vertex_coord);
                            subdomains[subDomain1].addVertex(v_bif);

                            counter_new_vessels += 1;
                        }
                    }
                }
            }
        }
    }

    if (counter_new_vessels > 0) {
        grid_copy->preGrow();

        // check mightVanish
        for (const auto& element : elements(grid_copy->levelGridView(0))) checkHierarchy(element);

        // Actually grow the network
        bool newElementGenerated = grid_copy->grow();
        if (!newElementGenerated) DUNE_THROW(Dune::Exception, "grid.grow() did not insert any new element");

        // Close growth process
        grid_copy->postGrow();
        auto gridView_copy = grid_copy->leafGridView();

        // Update parameter vectors
        G.resize(gridView_copy.size(0), 0.0);
        R.resize(gridView_copy.size(0), 0.0);
        mu.resize(gridView_copy.size(0), 0.0);
        u.resize(gridView_copy.size(1), 0.0);
        pO2[_0].resize(gridView_copy.size(1), 0.0);

        for (int i = 0; i < R.N(); i++) {
            R[i] = R_supp[i];
            G[i] = G_supp[i];
            mu[i] = mu_supp[i];
        }

        for (int i = 0; i < u.N(); i++) {
            u[i] = u_supp[i];
            pO2[_0][i] = pO2_supp[i];
        }

        return gridView_copy;
    }

    return gv1d;
}
