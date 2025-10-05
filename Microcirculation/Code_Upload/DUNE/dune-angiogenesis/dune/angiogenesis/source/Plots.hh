template <typename GV1D, typename Vector>
void plotNetwork(const GV1D& gv1d, Vector& R, Vector& G, Vector& mu, Vector& P, const std::string path,
                 std::string networkType, int day = 0, int refinement_level_1d = 0) {
    Vector u_mmHg = P;
    u_mmHg /= 133.32239;

    Dune::VTKWriter<GV1D> vtkWriter1d(gv1d, Dune::VTK::nonconforming);
    vtkWriter1d.addCellData(R, "R");
    vtkWriter1d.addCellData(G, "G");
    vtkWriter1d.addCellData(mu, "mu");
    vtkWriter1d.addVertexData(u_mmHg, "Pressure [mmHg]");

    std::string name = path;
    name.append(networkType);
    name.append("-day-");
    name += std::to_string(day);
    vtkWriter1d.write(name);
}

template <typename GV, typename Vector>
void plot1DSolution(const GV& gv, Vector& u, Vector& R, const std::string path, const std::string type, int day = 0,
                    int timestep = 0, bool exactSolution = false) {
    typedef double RF;
    typedef typename GV::template Codim<0>::Iterator LeafIterator;
    typedef typename GV::IntersectionIterator IntersectionIterator;
    typedef typename GV::IndexSet LeafIndexSet;
    const LeafIndexSet& set = gv.indexSet();
    typedef typename GV::Grid GridType;

    const int dimnetwork = GV::dimension;
    int numberOfElements = gv.size(0);

    Dune::VTKWriter<GV> vtkWriter(gv, Dune::VTK::nonconforming);

    Vector Rmu = R;
    Rmu *= 1E6;
    vtkWriter.addCellData(Rmu, "R");

    Vector u_mmHg = u;
    if (!exactSolution) u_mmHg /= 133.32239;

    if (type == "flow") {
        vtkWriter.addVertexData(u_mmHg, "Pressure [mmHg]");
    } else if (type == "flow_explicit") {
        vtkWriter.addVertexData(u_mmHg, "Pressure [mmHg]");
    } else if (type == "Glu" || type == "Glu_explicit"){
        vtkWriter.addVertexData(u_mmHg, "Glucose [mM]");
    } else if (type == "PUFA" || type == "PUFA_explicit"){
        vtkWriter.addVertexData(u_mmHg, "PUFA [mM]");
    }
    else {
        vtkWriter.addVertexData(u_mmHg, "PO_2 [mmHg]");
    }

    std::string name = path;
    name.append(type);
    name.append("-1D");
    if (type == "flow") {
        name.append("-day-");
        name += std::to_string(day);
    } else if (type == "flow_explicit") {
        name.append("-step-");
        name += std::to_string(day);
    } else {
        name.append("-day-");
        name += std::to_string(day);
        name.append("-timestep-");
        name += std::to_string(timestep);
    }
    vtkWriter.write(name);
}

template <typename GV, typename Vector>
void plot3DSolution(const GV& gv, Vector& u, const std::string path, const std::string type, int day = 0,
                    int timestep = 0, bool exactSolution = false) {
    Dune::VTKWriter<GV> vtkWriter(gv, Dune::VTK::conforming);

    Vector u_mmHg = u;
    if (!exactSolution) u_mmHg /= 133.32239;

    if (type == "flow") {
        vtkWriter.addCellData(u_mmHg, "Pressure [mmHg]");
    } else if (type == "Glu" || type == "Glu_explicit"){
        vtkWriter.addCellData(u_mmHg, "Glucose [mM]");
    } else if (type == "PUFA" || type == "PUFA_explicit"){
        vtkWriter.addCellData(u_mmHg, "PUFA [mM]");
    } else if (type == "TAF") {
        vtkWriter.addCellData(u, "TAF [pM]");
    } else if (type == "flow_explicit") {
        vtkWriter.addCellData(u_mmHg, "Pressure [mmHg]");
    } else {
        vtkWriter.addCellData(u_mmHg, "PO_2 [mmHg]");
    }

    std::string name = path;
    name.append(type);
    name.append("-3D");
    if (type == "flow") {
        name.append("-day-");
        name += std::to_string(day);
    } else if (type == "flow_explicit") {
        name.append("-step-");
        name += std::to_string(day);
    } else if (type == "TAF") {
        name.append("-level-");
        int level_3d = gv.grid().maxLevel();
        name += std::to_string(level_3d);
        name.append("-day-");
        name += std::to_string(day);
    } else {
        name.append("-day-");
        name += std::to_string(day);
        name.append("-timestep-");
        name += std::to_string(timestep);
    }
    vtkWriter.write(name);
}

template <typename GV>
void plot3DMesh(const GV& gv, const std::string path, const std::string type) {
    Dune::VTKWriter<GV> vtkWriter(gv, Dune::VTK::conforming);
    std::string name = path;
    name.append(type);
    name.append("-level-");
    int level_3d = gv.grid().maxLevel();
    name += std::to_string(level_3d);
    vtkWriter.write(name);
}
