template <typename GV, typename Vector>
void reportDataNetwork(const GV& gv, const std::string path, Vector& R) {
    typedef double RF;

    std::ios_base::openmode mode = std::ios_base::app;
    std::ofstream out(path, mode);

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
