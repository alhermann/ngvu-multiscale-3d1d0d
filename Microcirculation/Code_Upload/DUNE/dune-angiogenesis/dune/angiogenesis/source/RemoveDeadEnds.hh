template <typename GV, typename Vector, typename RF>
auto removeDeadEnds(const GV& gv, Vector& radiiVector, Vector& permeabilityVector, Vector& nodePressureVector,
                    Dune::FieldVector<RF, 3> O, Dune::FieldVector<RF, 3> L) {
    typedef typename GV::Grid NetworkGrid;
    const int dimworld = GV::dimensionworld;
    const int dimnetwork = GV::dimension;

    Dune::GridFactory<NetworkGrid> factory;
    Dune::GeometryType type(1);

    for (const auto& vertex : vertices(gv)) factory.insertVertex(vertex.geometry().center());

    for (const auto& element : elements(gv)) {
        std::vector<unsigned int> cornerIDs(2);
        cornerIDs[0] = gv.indexSet().subIndex(element, 0, dimnetwork);
        cornerIDs[1] = gv.indexSet().subIndex(element, 1, dimnetwork);
        factory.insertElement(type, cornerIDs);
    }

    auto grid_copy = factory.createGrid();

    Dune::PersistentContainer<NetworkGrid, RF> permeabilities(*grid_copy, 0);
    Dune::PersistentContainer<NetworkGrid, RF> radii(*grid_copy, 0);
    Dune::PersistentContainer<NetworkGrid, RF> pressures(*grid_copy, 1);

    for (const auto& element : elements(grid_copy->leafGridView())) {
        auto eIdx = grid_copy->leafGridView().indexSet().index(element);
        permeabilities[element] = permeabilityVector[eIdx];
        radii[element] = radiiVector[eIdx];
    }

    for (const auto& vertex : vertices(grid_copy->leafGridView())) {
        pressures[vertex] = nodePressureVector[grid_copy->leafGridView().indexSet().index(vertex)];
    }

    std::vector<bool> elementHasToBeEliminated(grid_copy->size(0), false);
    for (const auto& element : elements(grid_copy->leafGridView())) {
        auto eIdx = grid_copy->leafGridView().indexSet().index(element);
        for (const auto& intersection : intersections(grid_copy->leafGridView(), element)) {
            if (intersection.boundary()) {
                Dune::FieldVector<RF, dimworld> vertex_coord = intersection.geometry().center();

                bool atBoundary = isPointAtBoundary(vertex_coord, O, L);

                if (!atBoundary) {
                    elementHasToBeEliminated[eIdx] = true;
                }
            }
        }
    }

    // Start removing single elements in the grid
    grid_copy->preGrow();

    for (const auto& element : elements(grid_copy->leafGridView())) {
        auto eIdx = grid_copy->leafGridView().indexSet().index(element);
        if (elementHasToBeEliminated[eIdx]) {
            grid_copy->removeElement(element);
        }
    }

    // Actually remove the elements from the grid
    grid_copy->grow();

    // Recover element values
    Vector newPermeabilities, newRadii;
    newPermeabilities.resize(grid_copy->leafGridView().size(0));
    newRadii.resize(grid_copy->leafGridView().size(0));
    for (const auto& element : elements(grid_copy->leafGridView())) {
        auto eIdx = grid_copy->leafGridView().indexSet().index(element);
        newPermeabilities[eIdx] = permeabilities[element];
        newRadii[eIdx] = radii[element];
    }

    // Recover pressures
    Vector newPressures;
    newPressures.resize(grid_copy->leafGridView().size(1));
    for (const auto& vertex : vertices(grid_copy->leafGridView())) {
        auto eIdx = grid_copy->leafGridView().indexSet().index(vertex);
        newPressures[eIdx] = pressures[vertex];
    }

    // Close reduction process
    grid_copy->postGrow();

    radiiVector = newRadii;
    permeabilityVector = newPermeabilities;
    nodePressureVector = newPressures;

    auto gridView_copy = grid_copy->leafGridView();

    return gridView_copy;
}
