/*
 * SupportFunctions.hh
 */

#pragma once

#include "boundingBox/boundingboxtree.hh"

bool fileExists(const std::string& filename) {
    struct stat buf;
    return (stat(filename.c_str(), &buf) == 0);
}

/**
 * For each vessel i=1..R_size:
 *   - Wait until <directory>/NVU_Vessel_i.csv exists
 *   - Keep polling until the file contains a row with time >= t (seconds)
 *   - Interpolate radius at t (unchanged linear interpolation)
 * Returns radii in radii_vec as [R_size x 1] (µm).
 */
template <typename Matrix>
void waitForFileAndInterpolate(const std::string& directory, double t, Matrix& radii_vec, std::size_t R_size)
{
    radii_vec.resize(R_size, 1);

    auto sleep_a_bit = []() {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    };
    const double eps_time = 1e-12;

    for (std::size_t vesselIndex = 1; vesselIndex <= R_size; ++vesselIndex) {
        // Build filename NVU_Vessel_X.csv
        std::ostringstream oss;
        oss << directory << "/NVU_Vessel_" << vesselIndex << ".csv";
        const std::string filename = oss.str();

        // 1) Wait for file to exist
        while (!fileExists(filename)) {
            sleep_a_bit();
        }

        double interpolatedRadius_um = 0.0;
        bool have_interpolated = false;

        // 2) Poll until there is data with time >= t; then interpolate
        while (!have_interpolated) {
            std::ifstream infile(filename);
            if (!infile.is_open()) {
                sleep_a_bit();
                continue;
            }

            std::vector<std::pair<double, double>> timeRadiusData;
            timeRadiusData.reserve(1024);

            std::string line;
            while (std::getline(infile, line)) {
                if (line.empty()) continue;
                std::stringstream ss(line);
                double time_val, radius_um;
                char comma;
                if ((ss >> time_val >> comma >> radius_um) && (comma == ',')
                    && std::isfinite(time_val) && std::isfinite(radius_um)) {
                    timeRadiusData.emplace_back(time_val, radius_um);
                }
            }
            infile.close();

            if (timeRadiusData.empty()) {
                sleep_a_bit();
                continue;
            }

            std::sort(timeRadiusData.begin(), timeRadiusData.end(),
                      [](const auto& a, const auto& b) { return a.first < b.first; });

            // Not enough yet? keep waiting
            if (timeRadiusData.back().first + eps_time < t) {
                sleep_a_bit();
                continue;
            }

            // === Interpolation (unchanged) ===
            auto it = std::lower_bound(
                timeRadiusData.begin(), timeRadiusData.end(),
                std::make_pair(t, 0.0),
                [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
                    return a.first < b.first;
                });

            if (it == timeRadiusData.begin()) {
                interpolatedRadius_um = it->second;
                have_interpolated = true;
            } else if (it == timeRadiusData.end()) {
                interpolatedRadius_um = timeRadiusData.back().second;
                have_interpolated = true;
            } else {
                const auto [t2, r2] = *it;
                const auto [t1, r1] = *(it - 1);
                const double denom = (t2 - t1);
                interpolatedRadius_um = (std::abs(denom) < eps_time) ? r2
                    : r1 + (r2 - r1) * ((t - t1) / denom);
                have_interpolated = true;
            }
            // === end interpolation ===
        }

        radii_vec[vesselIndex - 1][0] = interpolatedRadius_um; // µm
    }
}

template <typename Vector>
void waitForFileAndRead(const std::string& filename, Vector &radii_vec, size_t R_size) {
    const std::string lockfile = filename + ".lock";

    // Check if the desired file exists before attempting to read
    while (!fileExists(filename)) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }

    // Poll until lock file no longer exists
    while (fileExists(lockfile)) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }

    // Check if the desired file exists before attempting to read
    if (!fileExists(filename)) {
        std::cerr << "File does not exist after lock removal: " << filename << std::endl;
        return;
    }

    // Proceed to read the main file once lock is released
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error opening the file for reading." << std::endl;
        return;
    }

    // Test print the contents of the file
    // std::string line;
    // while (std::getline(infile, line)) {
    //     std::cout << line << std::endl;
    // }

    int rowCount = 0;
    std::string line;
    while (std::getline(infile, line)){
        ++rowCount;
    }

    // Reset file pointer to the beginning
    infile.clear();  // Clear EOF flag
    infile.seekg(0, std::ios::beg);

    // Store the outputs from the neurons file in a dynamic DUNE matrix
    radii_vec.resize(rowCount, R_size);
    int i = 0;
    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::string token;
        int j = 0;
        while (std::getline(ss, token, ',')) {
            radii_vec[i][j] = std::stod(token);
            ++j;
        }
        ++i;
    }
    infile.close();
}

template <typename RF>
void loadFlowData(const std::string& filename, std::vector<Dune::FieldVector<RF, 1>>& flow) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << "\n";
        return;
    }

    std::string line;
    bool isFirstLine = true;
    RF unit_fac = 1e-6; // convert units from cm³/s to m³/s            

    // Read the file line by line
    while (std::getline(file, line)) {
        if (isFirstLine) {
            // Skip the header line
            isFirstLine = false;
            continue;
        }

        std::istringstream ss(line);
        std::string token;
        int columnIndex = 0;
        RF flowValue = 0.0;

        // Parse each column in the line
        while (std::getline(ss, token, ',')) {
            if (columnIndex == 2) { // 2 corresponds to "Flow q" column
                flowValue = std::stod(token); // Convert to RF (double by default)
                break;
            }
            ++columnIndex;
        }

        // Add the value to the flow vector
        flow.emplace_back(Dune::FieldVector<RF, 1>{flowValue * unit_fac}); 
    }

    file.close();
    std::cout << "Loaded " << flow.size() << " flow data points from " << filename << ".\n";
}

template <typename MatrixDType>
void writeMatrixToFile(MatrixDType& matrix, const std::string& filename) {
    // Create the lock file name
    std::string lockFileName = filename + ".lock";

    // Create the lock file
    {
        std::ofstream lockFile(lockFileName);
        if (!lockFile) {
            std::cerr << "Error creating lock file: " << lockFileName << std::endl;
            return;
        }
        // Close the lock file immediately after creation
        lockFile.close();
    }
    
    // Open an output file stream
    std::ofstream outFile(filename);

    // Check if the file was opened successfully
    if (!outFile) {
        std::cerr << "Error: Unable to open file for writing: " << filename << std::endl;
        return;  // or handle the error as appropriate
    }

    // Iterate over the matrix and write each entry to the file
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            outFile << matrix[i][j];

            // Add a comma after each entry, except the last one in the row
            if (j != matrix[i].size() - 1) {
                outFile << ",";
            }
        }
        outFile << "\n";  // New line at the end of each row
    }

    // Close the file
    outFile.close();

    // Remove the lock file
    std::remove(lockFileName.c_str());
}

template <typename Element>
Element getCoarseElement(Element e) {
    for (auto levelIdx = e.level(); levelIdx != 0; --levelIdx) e = e.father();
    return e;
}

template <typename GV>
std::pair<Dune::FieldVector<double, 3>, Dune::FieldVector<double, 3> > determineExtrema(const GV& gv,
                                                                                        bool enlargeDomain = false) {
    typedef double RF;
    const int dimworld = 3;
    typedef typename GV::template Codim<0>::Iterator LeafIterator;
    Dune::FieldVector<double, 3> A, B;
    const LeafIterator itend = gv.template end<0>();

    // Set first definition of vertices A and B
    for (LeafIterator it = gv.template begin<0>(); it != itend; ++it) {
        A = it->geometry().corner(0);
        B = it->geometry().corner(1);
        break;
    }

    // Find the minimal and maximal vertices
    for (LeafIterator it = gv.template begin<0>(); it != itend; ++it) {
        Dune::FieldVector<RF, 3> vertex0 = it->geometry().corner(0), vertex1 = it->geometry().corner(1);

        // Test if some component of vertex0 are minimal
        if (vertex0[0] < A[0]) A[0] = vertex0[0];
        if (vertex0[1] < A[1]) A[1] = vertex0[1];
        if (vertex0[2] < A[2]) A[2] = vertex0[2];

        // Test if some component of vertex1 are minimal
        if (vertex1[0] < A[0]) A[0] = vertex1[0];
        if (vertex1[1] < A[1]) A[1] = vertex1[1];
        if (vertex1[2] < A[2]) A[2] = vertex1[2];

        // Test if some component of vertex0 are maximal
        if (vertex0[0] > B[0]) B[0] = vertex0[0];
        if (vertex0[1] > B[1]) B[1] = vertex0[1];
        if (vertex0[2] > B[2]) B[2] = vertex0[2];

        // Test if some component of vertex0 are maximal
        if (vertex1[0] > B[0]) B[0] = vertex1[0];
        if (vertex1[1] > B[1]) B[1] = vertex1[1];
        if (vertex1[2] > B[2]) B[2] = vertex1[2];
    }

    // Enlarge domain of 20% in each direction
    if (enlargeDomain) {
        Dune::FieldVector<RF, dimworld> diff = B - A;
        diff /= 10.0;
        B += diff;
        A -= diff;
    }

    return {A, B};
}

template <typename BlockMatrix, typename Matrix>
void packBlockMatrix(BlockMatrix& A, Matrix& B) {
    using namespace Dune::Indices;
    using namespace Dune::Hybrid;

    typedef typename Matrix::RowIterator RowIterator;
    typedef typename Matrix::ColIterator ColIterator;

    std::vector<std::set<int> > pattern;

    Dune::FieldVector<int, Dune::Hybrid::size(A) + 1> N(0);
    int dofs = 0;
    forEach(integralRange(Dune::Hybrid::size(A)), [&](const auto i) {
        dofs += A[i][i].N();
        N[i + 1] = A[i][i].N();
    });
    for (int i = 1; i < Dune::Hybrid::size(A) + 1; i++) {
        N[i] += N[i - 1];
    }

    B.setSize(dofs, dofs);
    pattern.resize(dofs);

    forEach(integralRange(Dune::Hybrid::size(A)), [&](const auto i) {
        forEach(integralRange(Dune::Hybrid::size(A[i])), [&](const auto j) {
            for (RowIterator row = A[i][j].begin(); row != A[i][j].end(); ++row)
                for (ColIterator col = row->begin(); col != row->end(); ++col) {
                    pattern[N[i] + row.index()].insert(N[j] + col.index());
                }
        });
    });

    B.setBuildMode(Matrix::random);
    for (int i = 0; i < dofs; i++) {
        B.setrowsize(i, pattern[i].size());
    }
    B.endrowsizes();
    for (int i = 0; i < dofs; i++) {
        std::template set<int>::iterator setend = pattern[i].end();
        for (std::template set<int>::iterator setit = pattern[i].begin(); setit != setend; ++setit) {
            B.addindex(i, *setit);
        }
    }
    B.endindices();
    B = 0.0;

    forEach(integralRange(Dune::Hybrid::size(A)), [&](const auto i) {
        forEach(integralRange(Dune::Hybrid::size(A[i])), [&](const auto j) {
            for (RowIterator row = A[i][j].begin(); row != A[i][j].end(); ++row)
                for (ColIterator col = row->begin(); col != row->end(); ++col) {
                    B[N[i] + row.index()][N[j] + col.index()] = A[i][j][row.index()][col.index()];
                }
        });
    });
}

template <typename BlockVector, typename Vector>
void packBlockVector(BlockVector& b, Vector& F) {
    using namespace Dune::Indices;
    using namespace Dune::Hybrid;

    // Pack all blocks together
    Dune::FieldVector<int, Dune::Hybrid::size(b) + 1> N(0);
    int dofs = 0;
    forEach(integralRange(Dune::Hybrid::size(b)), [&](const auto i) {
        dofs += b[i].N();
        N[i + 1] = b[i].N();
    });
    F.resize(dofs);
    F = 0.0;
    for (int i = 1; i < Dune::Hybrid::size(b) + 1; i++) {
        N[i] += N[i - 1];
    }

    forEach(integralRange(Dune::Hybrid::size(b)), [&](const auto i) {
        for (int j = 0; j < b[i].N(); j++) {
            if (b[i][j] != 0.0) {
                F[N[i] + j] = b[i][j];
            }
        }
    });
}

template <typename Vector, typename BlockVector>
void spackVector(Vector& U, BlockVector& u) {
    using namespace Dune::Indices;
    using namespace Dune::Hybrid;

    // Subdivide the blocks
    Dune::FieldVector<int, Dune::Hybrid::size(u) + 1> N(0);
    forEach(integralRange(Dune::Hybrid::size(u)), [&](const auto i) { N[i + 1] = u[i].N(); });
    for (int i = 1; i < Dune::Hybrid::size(u) + 1; i++) {
        N[i] += N[i - 1];
    }

    forEach(integralRange(Dune::Hybrid::size(u)), [&](const auto i) {
        for (int j = 0; j < u[i].N(); j++) {
            u[i][j] = U[N[i] + j];
        }
    });
}


template <typename RF>
RF computeInVivoViscosity(RF r) {
    r *= 2.0 / 1e-6;
    RF mu = 1e-3;  // Initializazion as the plasma viscosity
    RF mu_h = 6.0 * std::exp(-0.085 * r) + 3.2 - 2.44 * std::exp(-0.06 * std::pow(r, 0.645));

    // RF H = 0.45; // Discarge hematocrit
    // RF C = ( 0.8+std::exp(-0.075*r) )*( -1.0+1.0/(1.0+std::pow(r,12.0)*1.0e-11) )
    // + 1.0/(1.0+std::pow(r,12.0)*1.0e-11);

    RF f = r / (r - 1.1);
    RF factor = std::pow(f, 2.0);
    mu *= factor * (1.0 + (mu_h - 1.0) * 1.0 * factor);
    return mu;
}

template <typename NetworkGV, typename RF, typename Vector>
void recomputeGv(NetworkGV& networkgv, const Vector& R, Vector& G, RF rho_bl) {
    // Assuming the size of R, G, and mu are the same and correspond to each vessel
    /*for (size_t i = 0; i < R.size(); ++i) {
        RF radius = R[i]; // Access the radius for the ith vessel
        RF viscosity = computeInVivoViscosity(radius);

        // Recompute conductance using the Hagen-Poiseuille equation
        RF conductance = (rho_bl * M_PI * std::pow(radius, 4)) / (8.0 * viscosity);

        // Update G with the new conductance value
        G[i] = conductance;
    }*/

    // Iterate over all elements in the grid view
    for (const auto& element : elements(networkgv)) {
        auto idxElement = networkgv.indexSet().index(element); // Get index of the current element

        RF radius = R[idxElement]; // Access the radius for the current element
        RF viscosity = computeInVivoViscosity(radius); // Access the viscosity for the current element
        RF length = element.geometry().volume(); // In 1D, volume() gives the length of the element

        // Recompute conductance using the updated radius and divide by the element's length
        RF conductance = (rho_bl * M_PI * std::pow(radius, 4)) / (8.0 * viscosity * length);

        // Update G with the new conductance value
        G[idxElement] = conductance;
    }
}

template <typename RF, int dim>
Dune::FieldVector<RF, dim> findPointOnCircle(Dune::FieldVector<RF, dim> p1, Dune::FieldVector<RF, dim> p2,
                                             Dune::FieldVector<RF, dim> xc, RF radius, RF theta) {
    RF c = std::cos(2.0 * M_PI * theta), s = std::sin(2.0 * M_PI * theta);
    Dune::FieldVector<RF, dim> n = (p2 - p1), result(0.0), rotator(0.0);
    n /= (p2 - p1).two_norm();
    rotator[0] = n[2];
    rotator[2] = -n[0];
    Dune::FieldMatrix<RF, 3, 3> R(0.0);
    R[0][0] = n[0] * n[0] * (1.0 - c) + c;
    R[1][1] = n[1] * n[1] * (1.0 - c) + c;
    R[2][2] = n[2] * n[2] * (1.0 - c) + c;
    R[0][1] = n[0] * n[1] * (1.0 - c) - n[2] * s;
    R[0][2] = n[0] * n[2] * (1.0 - c) + n[1] * s;
    R[1][0] = n[0] * n[1] * (1.0 - c) + n[2] * s;
    R[2][0] = n[0] * n[2] * (1.0 - c) - n[1] * s;
    R[1][2] = n[1] * n[2] * (1.0 - c) - n[0] * s;
    R[2][1] = n[1] * n[2] * (1.0 - c) + n[0] * s;

    R.mv(rotator, result);
    result *= radius;
    result += xc;
    return result;
}

template <typename RF, int dim>
Dune::FieldVector<RF, dim> pointOnSegment(Dune::FieldVector<RF, dim> p1, Dune::FieldVector<RF, dim> p2,
                                          Dune::FieldVector<RF, 1> s) {
    Dune::FieldVector<RF, dim> x(0.0);
    x += p2 - p1;
    x *= s[0] / (p2 - p1).two_norm();
    x += p1;
    return x;
}

template <typename Vector>
bool alreadyPresent(Vector alreadyMet, int n, int& indexFound) {
    bool value = false;
    for (int i = 0; i < alreadyMet.size(); i++)
        if (n == alreadyMet[i]) {
            value = true;
            indexFound = i;
            break;
        }

    return value;
}

template <typename Vector>
bool equilibriumReached(Vector& pO2, Vector& pO2_old, Vector& Glu, Vector& Glu_old, Vector& Pufa, Vector& Pufa_old) {
    Vector pO2_ = pO2;
    pO2_ -= pO2_old;

    Vector Glu_ = Glu;
    Glu_ -= Glu_old;

    Vector Pufa_ = Pufa;
    Pufa_ -= Pufa_old;

    double relativeError = 0.0;
    for (int i = 0; i < pO2_.N(); i++) pO2_[i] /= pO2[i];
    double relativeErrorO2 = 0.0;
    for (int i = 0; i < pO2_.N(); i++) relativeErrorO2 += std::abs(pO2_[i]);
    relativeErrorO2 /= pO2_.N();

    for (int i = 0; i < Glu_.N(); i++) Glu_[i] /= Glu[i];
    double relativeErrorGlu = 0.0;
    for (int i = 0; i < Glu_.N(); i++) relativeErrorGlu += std::abs(Glu_[i]);
    relativeErrorGlu /= Glu_.N();

    for (int i = 0; i < Pufa_.N(); i++) Pufa_[i] /= Pufa[i];
    double relativeErrorPufa = 0.0;
    for (int i = 0; i < Pufa_.N(); i++) relativeErrorPufa += std::abs(Pufa_[i]);
    relativeErrorPufa /= Pufa_.N();

    relativeError = std::max(std::max(relativeErrorGlu, relativeErrorO2), relativeErrorPufa);

    std::cout << "Relative difference between old and new solution = " << relativeError * 100 << " %"
              << "\n";
    if (relativeError < 5.0e-2) {
        return true;
    } else {
        return false;
    }
}

template <typename Vector>
void plot1DSolutionToCSV(Vector& u, Vector& R, const std::string& path, int step) {
    std::vector<double> u_cell(R.size()); // Adjust size if needed

    // Convert u from vertex to cell data
    // This is a placeholder - actual conversion logic depends on your grid and data structure
    for (int i = 0; i < u.size(); ++i) {
        // Simple average conversion; replace with your actual conversion logic
        u_cell[i] = (u[i] + u[i + 1]) / 2.0;
    }
    // Append last value or perform a suitable operation for the last cell if necessary
    u_cell[u.size()] = (u[u.size()] + u[u.size() - 1]) / 2.0; // This is a simplification

    // Adjust path and file name as necessary
    std::string filePath = path + "-step-" + std::to_string(step) + ".csv";
    
    // Create the lock file name
    std::string lockFileName = filePath + ".lock";

    // Create the lock file
    {
        std::ofstream lockFile(lockFileName);
        if (!lockFile) {
            std::cerr << "Error creating lock file: " << lockFileName << std::endl;
            return;
        }
        // Close the lock file immediately after creation
        lockFile.close();
    }

    std::ofstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Failed to open file for writing.\n";
        return;
    }
    
    // Write headers
    // file << "R,Pressure [mmHg]\n";
    
    // Write data
    for (size_t i = 0; i < R.size(); ++i) {
        file << R[i] << "," << u_cell[i] << "\n";
    }
    
    file.close();

    // Remove the lock file
    std::remove(lockFileName.c_str());
}


template <typename Vector>
void writeSolutionToFile(Vector v, std::string filename, std::string name, int step) {
    // Construct the filename
    filename.append(name);
    filename.append("-step-");
    filename += std::to_string(step);
    filename.append(".csv");

    // Create the lock file name
    std::string lockFileName = filename + ".lock";

    // Create the lock file
    {
        std::ofstream lockFile(lockFileName);
        if (!lockFile) {
            std::cerr << "Error creating lock file: " << lockFileName << std::endl;
            return;
        }
        // Close the lock file immediately after creation
        lockFile.close();
    }

    // Write the actual data
    std::ios_base::openmode mode = std::ios_base::trunc;
    std::ofstream out(filename, mode);
    for (int i = 0; i < v.N(); i++) {
        out << v[i] << " \n";
    }
    out.close();

    // Remove the lock file
    std::remove(lockFileName.c_str());
}

template <typename Vector>
void writeSolutionToFile(Vector v, std::string name, int step) {
    std::string filename = "result_vectors/";
    filename.append(name);
    filename.append("-step");
    filename += std::to_string(step);
    filename.append(".m");
    std::ios_base::openmode mode = std::ios_base::trunc;
    std::ofstream out(filename, mode);
    for (int i = 0; i < v.N(); i++) {
        out << i << " " << v[i] << " \n";
    }
    out.close();
}

template <typename Vector>
void readSolutionFromFile(Vector& v, int N, std::string name, int step) {
    // v.resize(N);
    std::string filename("result_vectors/");
    filename.append(name);
    filename.append("-step");
    filename += std::to_string(step);
    filename.append(".m");
    std::ostringstream out_filename;
    out_filename << filename.c_str();
    std::string file;
    file = out_filename.str().c_str();
    std::ifstream in_filename;
    in_filename.open(file.c_str());

    if (N == -1)
    {
        // Determine the number of lines in the file
        N = std::count(std::istreambuf_iterator<char>(in_filename), 
                        std::istreambuf_iterator<char>(), '\n');
        in_filename.clear(); // Clear EOF flag
        in_filename.seekg(0, std::ios::beg); // Move file pointer back to the beginning
    }
    
    v.resize(N);

    for (int i = 0; i < N; i++) {
        Dune::FieldVector<double, 2> supp(0.0);
        in_filename >> supp[0] >> supp[1];
        v[i] = supp[1];
    }

    in_filename.close();
}

template <typename GV1D, typename Vector>
void writeDGFFile(const GV1D& gv1d, Vector& R, Vector& p, int step, const std::string name) {
    std::string filename;
    filename.append(name);
    filename.append("-step-");
    filename += std::to_string(step);
    filename.append(".dgf");
    std::ios_base::openmode mode = std::ios_base::trunc;
    std::ofstream out(filename, mode);
    out << "DGF"
        << " \n";
    out << "Vertex"
        << " \n";
    out << "parameters 1"
        << " \n";
    for (const auto& vertex : vertices(gv1d)) {
        out << vertex.geometry().center() << " " << p[gv1d.indexSet().index(vertex)] << " \n";
    }
    out << "#"
        << " \n";
    out << "SIMPLEX"
        << " \n";
    out << "parameters 1"
        << " \n";
    for (const auto& element : elements(gv1d)) {
        out << gv1d.indexSet().subIndex(element, 0, 1) << " " << gv1d.indexSet().subIndex(element, 1, 1) << " "
            << R[gv1d.indexSet().index(element)] << "\n";
    }
    out << "#"
        << " \n";
    out << "BOUNDARYDOMAIN"
        << " \n";
    out << "default 1"
        << " \n";
    out << "#"
        << " \n";
    out.close();
}

template <typename GV3D, typename Vector>
double calculateAveragePressure(const GV3D& gv3d, Vector& p) {
    const int dimworld = GV3D::dimension;
    double value = 0.0;
    int n = 0;

    for (const auto& element : elements(gv3d)) {
        Dune::FieldVector<double, dimworld> x_c = element.geometry().center();
        value += p[gv3d.indexSet().index(element)];
        n += 1;
    }

    std::cout << "Average pressure in the tissue = " << (value / n) / 133.32239 << " [mmHg]\n";

    if (n == 0) {
        DUNE_THROW(Dune::Exception, "No elements in the FV-mesh are contained in the original element!");
    }

    return value / n;
}

template <typename RF>
int findSubomain(Dune::FieldVector<RF, 3> x, Dune::FieldVector<RF, 3> O, Dune::FieldVector<RF, 3> L, int NSubdomain) {
    int subDomain = 0;

    for (int j = 0; j < 3; j++) {
        RF diff = x[j] - O[j];
        RF edge_length = L[j] / NSubdomain;
        int number = 0;
        int factor = 1;

        for (int k = 0; k < j; k++) {
            factor *= NSubdomain;
        }

        for (int l = 0; l < NSubdomain; l++) {
            if (diff < (l + 1) * edge_length && l * edge_length < diff) {
                number = l;
                break;
            } else {
                number = 0;
            }
        }

        subDomain += number * factor;
    }

    return subDomain;
}

template <typename GV3D, typename Vector, typename RF>
std::vector<double> calculateLocalAveragePO2(const GV3D& gv3d, Dune::FieldVector<RF, 3> O, Dune::FieldVector<RF, 3> L,
                                             Vector& pO2, int NSubdomain) {
    const int dimworld = GV3D::dimension;
    std::vector<RF> values;
    std::vector<int> counter;

    for (int i = 0; i < NSubdomain * NSubdomain * NSubdomain; i++) {
        values.push_back(0.0);
        counter.push_back(0);
    }

    for (const auto& element : elements(gv3d)) {
        Dune::FieldVector<RF, dimworld> x_c = element.geometry().center();

        int subDomain = findSubomain(x_c, O, L, NSubdomain);

        values[subDomain] += pO2[gv3d.indexSet().index(element)];
        counter[subDomain] += 1;
    }

    for (int i = 0; i < NSubdomain * NSubdomain * NSubdomain; i++) {
        if (counter[i] != 0) {
            values[i] /= counter[i];
        } else {
            values[i] = 0;
        }

        std::cout << "Average PO2 pressure in the tissue = " << values[i] / 133.32239 << " [mmHg] "
                  << " Control volume i: " << i << "\n";
    }

    return values;
}

template <typename RF>
bool isPointAtBoundary(Dune::FieldVector<RF, 3> point, Dune::FieldVector<RF, 3> O, Dune::FieldVector<RF, 3> L) {
    bool atBoundary = false;

    for (int i = 0; i < 3; i++) {
        if (point[i] < O[i] + 1.0e-5 || point[i] > O[i] + L[i] - 1.0e-5) {
            atBoundary = true;
            break;
        }
    }

    return atBoundary;
}

template <typename GV1D, typename Vector, typename RF>
int countNumberOfLargeTerminalVessels(const GV1D& gv1d, Vector& R, Dune::FieldVector<RF, 3> O,
                                      Dune::FieldVector<RF, 3> L, RF threshold_LV) {
    typedef typename GV1D::Grid Grid;

    const int dimnetwork = GV1D::dimension;
    const int dimworld = 3;

    int numberOfLargeTerminalVessels = 0;

    // Initialize new network
    Dune::GridFactory<Grid> factory;
    Dune::GeometryType type(1);

    // Generate new network
    for (const auto& element : elements(gv1d)) {
        for (const auto& intersection : intersections(gv1d, element)) {
            if (intersection.boundary()) {
                int elementNumber = gv1d.indexSet().index(element);

                Dune::FieldVector<RF, dimworld> vertex_coord = intersection.geometry().center();

                bool atBoundary = isPointAtBoundary(vertex_coord, O, L);

                RF R_p = R[elementNumber];

                if (R_p > threshold_LV && !atBoundary) {
                    numberOfLargeTerminalVessels++;
                }
            }
        }
    }

    return numberOfLargeTerminalVessels;
}

template <typename GV1D, typename Vector, typename RF>
RF calculateOutletArea(const GV1D &gv1d, Vector &R, Vector &R0, Vector &P_init, RF mean_pressure_value){
    const int dimnetwork = GV1D::dimension;
    RF outletArea = 0.0;

    for (const auto& element1d : elements(gv1d)) {
        int index0 = gv1d.indexSet().subIndex(element1d, 0, dimnetwork),
            index1 = gv1d.indexSet().subIndex(element1d, 1, dimnetwork),
            elementNumber = gv1d.indexSet().index(element1d);

        for (const auto& intersection : intersections(gv1d, element1d)) {
            if (intersection.boundary()) {
                int index = gv1d.indexSet().subIndex(element1d, intersection.indexInInside(), dimnetwork);

                // Outlet boundary points less than pressure threshold
                if (R0[elementNumber] < 2.9e-6) {
                    outletArea += R[elementNumber] * R[elementNumber] * M_PI; 
                }
            }
        }
    }
    std::cout << "Total area: " << outletArea << std::endl; 
    return outletArea;
}

template <typename GV1D, typename Vector, typename RF>
int countNumberOfTerminalVessels(const GV1D& gv1d, Vector& R, Dune::FieldVector<RF, 3> O, Dune::FieldVector<RF, 3> L) {
    typedef typename GV1D::Grid Grid;

    const int dimnetwork = GV1D::dimension;
    const int dimworld = 3;

    int numberOfLargeTerminalVessels = 0;

    // Initialize new network
    Dune::GridFactory<Grid> factory;
    Dune::GeometryType type(1);

    // Generate new network
    for (const auto& element : elements(gv1d)) {
        for (const auto& intersection : intersections(gv1d, element)) {
            if (intersection.boundary()) {
                int elementNumber = gv1d.indexSet().index(element);

                Dune::FieldVector<RF, dimworld> vertex_coord = intersection.geometry().center();

                bool atBoundary = isPointAtBoundary(vertex_coord, O, L);

                RF R_p = R[elementNumber];

                if (!atBoundary) {
                    numberOfLargeTerminalVessels++;
                }
            }
        }
    }

    return numberOfLargeTerminalVessels;
}

template <typename Subdomain, typename RF>
bool testCollision(std::vector<Subdomain>& subdomains, int subDomain, Dune::FieldVector<RF, 3> p1,
                   Dune::FieldVector<RF, 3> d1, RF R_vessel) {
    bool isColliding = false;

    std::vector<std::vector<unsigned int> > indices = subdomains[subDomain].getIndices();
    std::vector<RF> radii = subdomains[subDomain].getRadii();
    std::vector<Dune::FieldVector<double, 3> > vertices = subdomains[subDomain].getVertices();

    int numberOfEdges = indices.size();

    for (int i = 0; i < numberOfEdges; i++) {
        int index1 = indices[i][0];
        int index2 = indices[i][1];

        RF R_other_vessel = radii[i];

        Dune::FieldVector<RF, 3> p2, d2;

        p2 = vertices[index1];
        d2 = vertices[index2];

        d2 = d2 - p2;

        if ((p2 - p1).two_norm() < 1.0e-10 || (p2 + d2 - p1).two_norm() < 1.0e-10 ||
            ((p2 + d2) - (p1 + d1)).two_norm() < 1.0e-10 || (p2 - (p1 + d1)).two_norm() < 1.0e-10) {
            continue;
        }

        Dune::FieldVector<RF, 3> point_l1, point_l2, delta_1, delta_2;

        for (int k = 1; k < 21; k++) {
            delta_1[0] = ((RF)k / 22.0 * d1[0]);
            delta_1[1] = ((RF)k / 22.0 * d1[1]);
            delta_1[2] = ((RF)k / 22.0 * d1[2]);

            for (int j = 1; j < 21; j++) {
                delta_2[0] = ((RF)j / 22.0 * d2[0]);
                delta_2[1] = ((RF)j / 22.0 * d2[1]);
                delta_2[2] = ((RF)j / 22.0 * d2[2]);

                point_l1 = p1 + delta_1;
                point_l2 = p2 + delta_2;

                RF dist = (point_l2 - point_l1).two_norm();

                if (dist < R_vessel + R_other_vessel) {
                    isColliding = true;
                    break;
                }

                if (isColliding == true) {
                    break;
                }
            }

            if (isColliding == true) {
                break;
            }
        }
    }

    return isColliding;
}

template <typename RF, int dimworld>
std::pair<bool, int> isVertexAlreadyPresent(std::vector<Dune::FieldVector<RF, dimworld> > matrixOfNodes,
                                            Dune::FieldVector<RF, dimworld> node) {
    bool found = false;
    double tol = 1e-12;
    int line = -1;
    for (int i = 0; i < matrixOfNodes.size(); i++) {
        if (std::abs(matrixOfNodes[i][0] - node[0]) < tol && std::abs(matrixOfNodes[i][1] - node[1]) < tol &&
            std::abs(matrixOfNodes[i][2] - node[2]) < tol) {
            found = true;
            line = i;
            break;
        }
    }
    return std::pair<bool, int>(found, line);
}

template <typename Vector>
bool isOutside(Vector v_bif, Vector O, Vector L) {
    bool value = false;
    if (v_bif[0] > L[0] || v_bif[1] > L[1] || v_bif[2] > L[2] || v_bif[0] < O[0] || v_bif[1] < O[1] ||
        v_bif[2] < O[2]) {
        value = true;
    }
    return value;
}

template <typename RF, int dim>
Dune::FieldVector<RF, dim> findPointOnCircleBifurcation(Dune::FieldVector<RF, dim> n, Dune::FieldVector<RF, dim> xc,
                                                        RF radius, RF phi, Dune::FieldVector<RF, dim> rotator) {
    Dune::FieldVector<RF, dim> result(0.0);
    n /= n.two_norm();
    rotator /= rotator.two_norm();

    RF c = std::cos(M_PI - phi), s = std::sin(M_PI - phi);
    Dune::FieldMatrix<RF, 3, 3> R(0.0);
    R[0][0] = n[0] * n[0] * (1.0 - c) + c;
    R[1][1] = n[1] * n[1] * (1.0 - c) + c;
    R[2][2] = n[2] * n[2] * (1.0 - c) + c;
    R[0][1] = n[0] * n[1] * (1.0 - c) - n[2] * s;
    R[0][2] = n[0] * n[2] * (1.0 - c) + n[1] * s;
    R[1][0] = n[0] * n[1] * (1.0 - c) + n[2] * s;
    R[2][0] = n[0] * n[2] * (1.0 - c) - n[1] * s;
    R[1][2] = n[1] * n[2] * (1.0 - c) - n[0] * s;
    R[2][1] = n[1] * n[2] * (1.0 - c) + n[0] * s;

    R.mv(rotator, result);
    result *= radius;
    result += xc;
    return result;
}

