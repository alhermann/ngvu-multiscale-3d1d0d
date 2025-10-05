class Subdomain {
   public:
    Subdomain() {}
    virtual ~Subdomain() {}

    double getRadius(int index) { return R[index]; }
    std::vector<double> getRadii() { return R; }
    void setRadius(int index, double radius) { R[index] = radius; }
    void addRadius(double radius) { R.push_back(radius); }

    Dune::FieldVector<double, 3> getVertex(int index) { return vertices[index]; }
    std::vector<Dune::FieldVector<double, 3> > getVertices() { return vertices; }
    void addVertex(Dune::FieldVector<double, 3> newVertex) { vertices.push_back(newVertex); }

    std::vector<std::vector<unsigned int> > getIndices() { return indices; }

    void addIndex() {
        int numberOfEdges = indices.size();
        std::vector<unsigned int> newIndices;
        newIndices.push_back(2 * numberOfEdges);
        newIndices.push_back(2 * numberOfEdges + 1);

        indices.push_back(newIndices);
    }

    std::vector<unsigned int> getIndex(int index) { return indices[index]; }

   private:
    std::vector<double> R;
    std::vector<Dune::FieldVector<double, 3> > vertices;
    std::vector<std::vector<unsigned int> > indices;
};
