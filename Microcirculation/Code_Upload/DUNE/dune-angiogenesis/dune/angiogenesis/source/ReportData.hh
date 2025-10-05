template <typename GV, typename Vector>
void reportDataNetwork(const GV& gv, const std::string path, Vector& R) {
    typedef double RF;

    std::string name = path;
    name.append("data_extracted_network.m");

    std::ios_base::openmode mode = std::ios_base::app;
    std::ofstream out(name, mode);

    for (const auto& element : elements(gv)) {
        auto eIdx = gv.indexSet().index(element);

        RF length = 0.0;

        for (const auto& intersection : intersections(gv, element)) {
            Dune::FieldVector<RF, 3> vertex_coord = intersection.geometry().center();
            Dune::FieldVector<RF, 3> vb = element.geometry().corner(1 - intersection.indexInInside());

            vb -= vertex_coord;
            length = vb.two_norm();
        }

        out << R[eIdx] << " " << length << " \n";
    }
}

void reportAveragePressue(const std::string path, double avg_PO2, double avg_Pressure) {
    std::string name = path;
    name.append("Average_PO2_Fluidpressure.txt");

    std::ios_base::openmode mode = std::ios_base::app;
    std::ofstream out(name, mode);

    out << avg_PO2 << " " << avg_Pressure << " \n";
}

void reportFluxIntoTissue(const std::string path, double flux) {
    std::string name = path;
    name.append("Flux_Into_Tissue.txt");

    std::ios_base::openmode mode = std::ios_base::app;
    std::ofstream out(name, mode);

    out << flux << " \n";
}
