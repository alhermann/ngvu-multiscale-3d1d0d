#pragma once

class SCV {
   public:
    SCV(int i) { identifier = i; }
    virtual ~SCV() {}

    int getSCVidentifier() { return identifier; }

    void saveRadius(double r_) { radius = r_; }
    double getRadius() { return radius; }

    void saveDetJac(double d_) { detjac = d_; }
    double getDetJac() { return detjac; }

    void saveWeight(double w_) { weights.push_back(w_); }
    std::vector<double> getWeights() { return weights; }

    void initializeStencil() {
        elements.push_back(std::vector<int>());
        lengths.push_back(std::vector<double>());
    }

    void saveStencilElement(int eIdx, double l_) {
        elements.back().push_back(eIdx);
        lengths.back().push_back(l_);
    }
    std::vector<int> getStencilElements(int i) { return elements[i]; }

    void addLength(double l_, int idx) { lengths.back()[idx] += l_; }
    std::vector<double> getLength(int i) { return lengths[i]; }

   private:
    int identifier = -1;
    std::vector<std::vector<int>> elements;
    std::vector<std::vector<double>> lengths;
    int centerElement;
    double radius = 0.0, detjac = 0.0;
    std::vector<double> weights;
};

class ControlVolume {
   public:
    ControlVolume() { containing3DElement = 0; }
    virtual ~ControlVolume() {}

    void addSCV(SCV* scv_) { scv.push_back(scv_); }
    std::vector<SCV*> getSubControlVolumes() { return scv; }

    bool checkSCV(int eIdx) {
        for (int k = 0; k < scv.size(); k++) {
            int identifier = scv[k]->getSCVidentifier();
            if (identifier == eIdx) {
                return true;
            }
        }
        return false;
    }

    void saveContaining3DElement(int eIdx) { containing3DElement = eIdx; }
    int getContaining3DElement() const { return containing3DElement; }

    void addVolume(double v) { volume += v; }
    double getVolume() const { return volume; }

    // Routines new coupling
    void save3DElementNewCoupling(int eIdx) { elements3DNewCoupling.push_back(eIdx); }
    std::vector<int> get3DElementsNewCoupling() const { return elements3DNewCoupling; }

    void saveNewWeight(double weight) { weights.push_back(weight); }
    void saveWeight(double weight, int idx) { weights[idx] += weight; }
    double getWeight(int idx) const { return weights[idx]; }
    std::vector<double> getWeights() const { return weights; }

   private:
    std::vector<SCV*> scv;
    int containing3DElement;
    double volume = 0.0;

    // Variables new coupling
    std::vector<double> weights;
    std::vector<int> elements3DNewCoupling;
};

template <typename GV1D, typename GV3D, typename Tree, typename Vector, typename RF>
std::vector<ControlVolume> determineStencils(const GV1D& gv1d, const GV3D& gv3d, Tree& tree, Vector& R, RF nQC,
                                             RF threshold = 7e-6) {
    const int dimworld = GV3D::dimension;
    const int dimnetwork = GV1D::dimension;

    std::vector<ControlVolume> controlVolumes(gv1d.size(dimnetwork));

    typedef typename GV3D::template Codim<0>::Entity Element;
    IndexToElementMap<GV3D> indexToElementMap(gv3d);
    indexToElementMap.resize(gv3d.size(0));
    for (auto&& element : elements(gv3d)) {
        unsigned int eIdx = gv3d.indexSet().index(element);
        indexToElementMap[eIdx] = element.seed();
    }

    // Create types for quadrature on subcontrol volumes
    Dune::GeometryType segment(Dune::GeometryType::simplex, 1);
    const Dune::QuadratureRule<RF, dimnetwork>& rule = Dune::QuadratureRules<RF, dimnetwork>::rule(segment, 4);

    // Determine coupled elements
    for (const auto& element1d : elements(gv1d)) {
        int eIdx1d = gv1d.indexSet().index(element1d);

        RF radius = R[eIdx1d];

        int index0 = gv1d.indexSet().subIndex(element1d, 0, dimnetwork);
        int index1 = gv1d.indexSet().subIndex(element1d, 1, dimnetwork);
        controlVolumes[index0].addVolume(0.5 * element1d.geometry().volume() * M_PI * radius * radius);
        controlVolumes[index1].addVolume(0.5 * element1d.geometry().volume() * M_PI * radius * radius);

        if (radius >= threshold) {
            continue;
        }

        Dune::FieldVector<RF, dimworld> x_m = element1d.geometry().center();
        RF surfaceArea = 2.0 * M_PI * radius;

        for (const auto& intersection1d : intersections(gv1d, element1d)) {
            int index1d = gv1d.indexSet().subIndex(element1d, intersection1d.indexInInside(), dimnetwork);
            Dune::FieldVector<RF, dimworld> x_i = intersection1d.geometry().center();
            Dune::FieldVector<RF, dimworld> direction = x_m;
            direction -= x_i;
            direction /= (x_m - x_i).two_norm();

            RF h = (x_m - x_i).two_norm() / nQC;

            for (int i = 0; i < nQC - 1; i++) {
                Dune::FieldVector<RF, dimworld> x_segment = x_i, update = direction;
                update *= (i + 1) * h;
                x_segment += update;

                RF H = h;
                if (i == 0 || i == nQC - 2) {
                    H = 1.5 * h;
                }

                for (int k = 0; k < nQC; k++) {
                    RF theta = k / nQC + 0.01;
                    Dune::FieldVector<RF, dimworld> pointCircle = findPointOnCircle(x_m, x_i, x_segment, radius, theta);
                    std::vector<unsigned int> elementsQuadOnCircle = tree.computeEntityCollisions(pointCircle);

                    if (elementsQuadOnCircle.size() == 0) {
                        continue;
                    }

                    Element elementQuadOnCircle = indexToElementMap.entity(elementsQuadOnCircle[0]);
                    int index3d = gv3d.indexSet().index(elementQuadOnCircle), indexFound = 0;

                    if (alreadyPresent(controlVolumes[index1d].get3DElementsNewCoupling(), index3d, indexFound)) {
                        controlVolumes[index1d].saveWeight(surfaceArea * H / nQC, indexFound);
                    } else {
                        controlVolumes[index1d].save3DElementNewCoupling(index3d);
                        controlVolumes[index1d].saveNewWeight(surfaceArea * H / nQC);
                    }
                }
            }
        }
    }

    return controlVolumes;
}
