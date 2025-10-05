template <typename GV1D, typename GV3D, typename CV, typename Vector, typename RF>
void computeFluxes(const GV1D& gv1d, const GV3D& gv3d, CV cvs, Vector& u1D, Vector& u3D, RF rho_int, RF Lp,
                   RF osmoticGradient, RF& flux_capillaries_tissue, RF& flux_tissue_capillaries) {
    const int dimnetwork = GV1D::dimension;
    const int dimworld = GV3D::dimension;

    for (auto it = cvs.begin(); it != cvs.end(); ++it) {
        // Get the control volume
        ControlVolume cv = *it;

        int index1d = it - cvs.begin();

        std::vector<int> couplingElements = cv.get3DElementsNewCoupling();
        int numberOfElements = couplingElements.size();
        std::vector<double> weights = cv.getWeights();

        RF pressure1d = u1D[index1d];

        // Flux constant
        RF factor = Lp * rho_int;

        for (int j = 0; j < numberOfElements; j++) {
            int index3d = couplingElements[j];
            RF pressure3d = u3D[index3d];

            if (pressure1d - pressure3d - osmoticGradient > 0.0) {
                flux_capillaries_tissue += factor * (pressure1d - pressure3d - osmoticGradient) * weights[j];
            } else {
                flux_tissue_capillaries += factor * (pressure1d - pressure3d - osmoticGradient) * weights[j];
            }

        }  // end loop over subcontrol volumes

    }  // end loop over control volumes
}