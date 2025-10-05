#include "ControlVolume.hh"

template <typename GV3D, typename GV1D, typename CV, typename BlockPattern>
void determinePattern(const GV1D& gv1d, const GV3D& gv3d, CV cvs, BlockPattern& blockPattern) {
    typedef double RF;
    using namespace Dune::Indices;

    const int dimnetwork = GV1D::dimension;
    const int dimworld = GV3D::dimension;
    int N1D = gv1d.size(dimnetwork), N3D = gv3d.size(0);

    blockPattern[_0][_0].clear();
    blockPattern[_0][_1].clear();
    blockPattern[_1][_0].clear();
    blockPattern[_1][_1].clear();

    blockPattern[_0][_0].resize(N1D);
    blockPattern[_0][_1].resize(N1D);
    blockPattern[_1][_0].resize(N3D);
    blockPattern[_1][_1].resize(N3D);

    // Determine standard pattern for the network
    int vertexsize = 2;
    for (const auto& element1d : elements(gv1d)) {
        for (int i = 0; i < vertexsize; i++) {
            int indexi = gv1d.indexSet().subIndex(element1d, i, dimnetwork);
            for (int j = 0; j < vertexsize; j++) {
                int indexj = gv1d.indexSet().subIndex(element1d, j, dimnetwork);
                blockPattern[_0][_0][indexi].insert(indexj);
            }
        }
    }

    std::cout << "Determination pattern for standard operator in network done. Computing pattern in tissue...\n";

    // Determine standard pattern for the tissue
    for (const auto& element : elements(gv3d)) {
        int indexi = gv3d.indexSet().index(element);
        blockPattern[_1][_1][indexi].insert(indexi);

        for (const auto& intersection : intersections(gv3d, element)) {
            if (!intersection.boundary()) {
                int indexj = gv3d.indexSet().index(intersection.outside());
                blockPattern[_1][_1][indexi].insert(indexj);
            }
        }
    }

    std::cout << "Determination pattern for standard operators done. Computing coupling...\n";

    // Determine pattern for the coupling
    for (auto it = cvs.begin(); it != cvs.end(); ++it) {
        // Get the control volume
        ControlVolume cv = *it;

        int index1d = it - cvs.begin();
        std::vector<int> couplingElements = cv.get3DElementsNewCoupling();
        int numberOfElements = couplingElements.size();

        for (int j = 0; j < numberOfElements; j++) {
            int index3d = couplingElements[j];
            blockPattern[_1][_0][index3d].insert(index1d);
            blockPattern[_0][_1][index1d].insert(index3d);
        }
    }
}
