#include "HarmonicPolyhedron.hpp"
#include "Tetrahedralize.hpp"
#include <iostream>

void loadPoly(const std::string& fname, std::vector<Eigen::Vector3d>& vertices, std::vector<std::vector<int>>& faces) {
    
    std::ifstream file(fname);
    
    int nv, nf;
    file >> nv;
    file >> nf;
    
    for(int i = 0; i < nv; ++i) {
        Eigen::Vector3d v;
        
        file >> v(0);
        file >> v(1);
        file >> v(2);
        
        vertices.push_back(v);
    }
    
    
    for(int i = 0; i < nf; ++i) {
        std::vector<int> face;
        int nfi;
        
        file >> nfi;
        
        for(int j = 0; j < nfi; ++j) {
            int id;
            file >> id;
            face.push_back(id);
        }
        
        faces.push_back(face);
    }
    
    file.close();
}

Edge::Edge(const int i_, const int j_) {
    if(i_ < j_) {
        i = i_;
        j = j_;
    } else {
        i = j_;
        j = i_;
    }
}

bool Edge::operator==(const Edge& e2) const {
    return i == e2.i && j == e2.j;
}

bool Edge::operator<(const Edge& e2) const {
    if(i < e2.i) return true;
    if(i > e2.i) return false;
    if(j < e2.j) return true;
    return false;
}


Normal::Normal()
: n(Eigen::Vector3d::Zero()) {}

void Normal::normalize() {
    n.normalize();
}

const Normal& Normal::operator+=(const Eigen::Vector3d& v) {
    n += v;
    return *this;
}
    

void HarmonicPolyhedron::dumpProbesAndKernels() {
    std::ofstream fk("kernels");
    for(const auto& v : kernelCenters) fk << ((v + mean) / scale).transpose() << std::endl;
    fk.close();
    
    std::ofstream fs("probes");
    
    // first face samples
    for(int i = 0; i < (int)harmonicFaces.size(); ++i) {
        auto samples3d = harmonicFaces[i].unproject(faceSamples2d[i]);
        
        for(auto& s : samples3d) {
            fs << ((s + mean) / scale).transpose() << std::endl;
        }
    }
    
    // second edge/vertex samples
    Eigen::MatrixXd vertexMatrix(vertices.size(), 3);
    for(int i = 0; i < (int)vertices.size(); ++i) vertexMatrix.row(i) = vertices[i].transpose();
    Eigen::MatrixXd evSamples = PEdgeSamples * vertexMatrix;
    
    for(int i = 0; i < evSamples.rows(); ++i) {
        fs << Eigen::RowVector3d((evSamples.row(i).transpose() + mean) / scale) << std::endl;
    }
    
    fs.close();
}

void HarmonicPolyhedron::dumpCrossection() {
    
    std::ofstream file("crossSamples");
    const int N = 1024;
    
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            Eigen::Vector3d s3d(-1. + 2. * i / (double)N, .0, -1. + 2. * j / (double)N);
            
            if(inside(s3d)) {
                file << evaluate(0, s3d) << " ";
            } else file << .0 << " ";
        }
        
        file << std::endl;
    }
    
    file.close();
}

bool HarmonicPolyhedron::insideConvex(const Eigen::Vector3d& p) {
    for(auto & harmonicFace : harmonicFaces) {
        if(p.dot(harmonicFace.getNormal()) > harmonicFace.getD()) return false;
    }

    return true;
}

bool HarmonicPolyhedron::inside(const Eigen::Vector3d& p) {
    // shoot ray in arbitrary direction (1, 0, 0) and intersect with all tris.t Add +/-1 depending on orientation
    
    Eigen::Vector3d dir;
    dir.setRandom();
    dir.normalize();
    
    int cnt = 0;
    
    for(auto& t : triPoints) {
        
        double dist;
        
        if(triangleRayIntersection(t[0], t[1], t[2], p, dir, dist)) {
            cnt += trianglePointOrientation(t[0], t[1], t[2], p) ? 1 : -1;
        }
    }
    
    return cnt != 0;
}

int HarmonicPolyhedron::globalToLocalIdMap(const int vid, const int fid) {
    for(int i = 0; i < (int)faces[fid].size(); ++i) {
        if(faces[fid][i] == vid) return i;
    }
    
    return -1;
}

Eigen::Vector3d HarmonicPolyhedron::evaluateGrad(const int node, const Eigen::Vector3d& p) {
    if(node < 0 || node >= coefficients.rows()) return Eigen::Vector3d::Zero();
    
    const Eigen::VectorXd c = coefficients.row(node);
    Eigen::Vector3d v = c.middleRows(nKernels, 3);
    
    for(int i = 0; i < nKernels; ++i) {
        const double nrm = (kernelCenters[i] - p).squaredNorm();
        v += c(i) / (sqrt(nrm) * nrm) * (kernelCenters[i] - p);
    }
    
    return v;
}

double HarmonicPolyhedron::evaluate(const int node, const Eigen::Vector3d& p) {
    if(node < 0 || node >= coefficients.rows()) return .0;
    
    const Eigen::VectorXd c = coefficients.row(node);
    double v = c.middleRows(nKernels, 3).dot(p) + c(nKernels + 3);
    
    for(int i = 0; i < nKernels; ++i) {
        v += c(i) / ((kernelCenters[i] - p).norm());
    }
    
    return v;
}

void HarmonicPolyhedron::initQuadrature() {
    
    std::vector<Eigen::MatrixXd> tets;
    tetrahedralize(vertices, faces, tets, 0.8);
   // std::cout << tets.size() << std::endl;
    tetQuad = TetrahedralQuadrature(tets);
}


void HarmonicPolyhedron::stiffnessMatrix(Eigen::MatrixXd& K) {
    
    const int n = (int)vertices.size();
    K.resize(n, n);
    
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j <= i; ++j) {
            double val = std::pow(1. / scale, 1) * tetQuad.apply([&](const Eigen::Vector3d p){return evaluateGrad(i, p).dot(evaluateGrad(j, p));});
            
            K(j, i) = K(i, j) = val;
        }
    }
}

void HarmonicPolyhedron::massMatrix(Eigen::MatrixXd& M) {
    
    const int n = (int)vertices.size();
    M.resize(n, n);
    
    for(int i = 0; i < n; ++i) {
          for(int j = 0; j <= i; ++j) {
              double val = std::pow(1. / scale, 3) * tetQuad.apply([&](const Eigen::Vector3d p){return evaluate(i, p) * evaluate(j, p);});
              
              M(j, i) = M(i, j) = val;
          }
      }
}


HarmonicPolyhedron::HarmonicPolyhedron() {}

HarmonicPolyhedron::HarmonicPolyhedron(const std::string& fname) {
    loadPoly(fname, vertices, faces);
    init();
}


HarmonicPolyhedron::HarmonicPolyhedron(const std::vector<Eigen::Vector3d>& vertices_, const std::vector<std::vector<int>>& faces_)
: vertices(vertices_), faces(faces_) {
    init();
}

void HarmonicPolyhedron::init() {
    
    using namespace std;
    
    ///////////////////////////////////////////////////
    // Center and rescale polyhedron
    
    mean.setZero();
    for(auto& v : vertices) mean += v;
    mean /= (double)vertices.size();
    scale = .0;
    
    for(auto& v : vertices) scale = std::max(scale, (v - mean).norm());
    scale = 1. / scale;
    
    for(auto& v : vertices) v = (v - mean) * scale;
    
    ///////////////////////////////////////////////////
    // Split faces into triangles: we use this for robust inside outside tests
    
    for(auto& f : faces) {
        const int n = (int)f.size();
        Eigen::Vector3d Mean(.0, .0, .0);
        for(int i = 0; i < n; ++i) {
            Mean += vertices[f[i]];
        }

        Mean /= (double)n;
        
        for(int i = 0; i < n; ++i) {
            triPoints.push_back(vector<Eigen::Vector3d>{vertices[f[i]], vertices[f[(i + 1) % n]] , Mean});
        }
    }
    
    ///////////////////////////////////////////////////
    // Sampling for probe points on faces and edges
    
    // sample faces for boundary conditions
    for(auto& f : faces) {
        Eigen::MatrixXd poly(f.size(), 3);
        for(int i = 0; i < (int)f.size(); ++i) poly.row(i) = vertices[f[i]];
        harmonicFaces.emplace_back(poly);
        faceSamples2d.push_back(harmonicFaces.back().generateRandomSamples(probesPerFace));
        nProbes += (int)faceSamples2d.back().size();
    }
    
    // generate unique edge and vertex samples.
    vector<Normal> vertexNormals(vertices.size());
    map<Edge, Normal> edgeNormals;
    
    for(int k = 0; k < (int)faces.size(); ++k) {
        auto& f = faces[k];
        const int n = (int)f.size();
        const Eigen::Vector3d areaVector = harmonicFaces[k].getScaledNormal();
        
        for(int i = 0; i < (int)f.size(); ++i) {
            vertexNormals[f[i]] += areaVector;
            edges.emplace_back(f[i], f[(i+1)%n]);
            edgeNormals[edges.back()] += areaVector;
        }
    }
    
    for(auto& n : edgeNormals) n.second.normalize();
    for(auto& n : vertexNormals) n.normalize();
    
    sort(edges.begin(), edges.end());
    edges.erase(unique(edges.begin(), edges.end()), edges.end());
    
    // Build sparse matrix mapping from vertex data to boundary conditions along edges/vertices.
    vector<Eigen::Triplet<double>> trip;
    int cnt = 0;
    
    for(int i = 0; i < (int)vertices.size(); ++i) {
        trip.emplace_back(cnt++, i, 1.);
    }
    
    for(auto & edge : edges) {
        for(int j = 1; j < probesPerEdge; ++j) {
            const double w = j / (double)probesPerEdge;
            trip.emplace_back(cnt, edge.i, 1. - w);
            trip.emplace_back(cnt, edge.j, w);
            ++cnt;
        }
    }
    
    PEdgeSamples.resize(cnt, (int)vertices.size());
    PEdgeSamples.setFromTriplets(trip.begin(), trip.end());
    nProbes += cnt;
    
    ///////////////////////////////////////////////////
    // Generate kernel centers by sampling all faces, edges and vertices
    
    for(auto& hf : harmonicFaces) {
        auto samples = hf.unproject(hf.generateRandomSamples(kernelsPerFace));
        for(auto& s : samples) kernelCenters.emplace_back(s + eps * hf.getNormal());
    }
    
    for(auto& e : edges) {
        const auto nrml = edgeNormals[e].n;
        
        for(int i = 1; i < kernelsPerEdge; ++i) {
            const double w = i / (double)kernelsPerEdge;
            kernelCenters.emplace_back((1. - w) * vertices[e.i] + w * vertices[e.j] + eps * nrml);
        }
    }
    
    for(int i = 0; i < (int)vertices.size(); ++i) {
        kernelCenters.emplace_back(vertices[i] + eps * vertexNormals[i].n);
    }
    
    
    nKernels = (int)kernelCenters.size();
    
    ///////////////////////////////////////////////////
    // Setup least squares system matrix
    
    Eigen::MatrixXd A(nProbes, nKernels + 4);
    cnt = 0;
    
    // first face samples
    for(int i = 0; i < (int)harmonicFaces.size(); ++i) {
        auto samples3d = harmonicFaces[i].unproject(faceSamples2d[i]);
        
        for(auto& s : samples3d) {
            for(int j = 0; j < (int)kernelCenters.size(); ++j) {
                A(cnt, j) = 1. / ((s - kernelCenters[j]).norm());
            }
            
            A.block(cnt, nKernels, 1, 3) = s.transpose();
            A(cnt, nKernels + 3) = 1.;
            
            ++cnt;
        }
    }
    
    // second edge/vertex samples
    Eigen::MatrixXd vertexMatrix(vertices.size(), 3);
    for(int i = 0; i < (int)vertices.size(); ++i) vertexMatrix.row(i) = vertices[i].transpose();
    Eigen::MatrixXd evSamples = PEdgeSamples * vertexMatrix;
    
    for(int i = 0; i < evSamples.rows(); ++i) {
        Eigen::Vector3d s = evSamples.row(i);
        
        for(int j = 0; j < (int)kernelCenters.size(); ++j) {
            A(cnt, j) = 1. / ((s - kernelCenters[j]).norm());
        }
        
        A.block(cnt, nKernels, 1, 3) = s.transpose();
        A(cnt, nKernels + 3) = 1.;
        
        ++cnt;
    }
    
    Eigen::LDLT<Eigen::MatrixXd> chol(A.transpose() * A);
    
    ////////////////////////////////////////////////////////////////////
    // Solve the system for each basis function to get the coefficients
    
    // for each vertex node build the correct rhs
    Eigen::VectorXd b(nProbes);
    coefficients.resize((int)vertices.size(), nKernels + 4);
    
    for(int i = 0; i < (int)vertices.size(); ++i) {
        cnt = 0;
        
        // first evaluate at face probes
        for(int j = 0; j < (int)harmonicFaces.size(); ++j) {
            const auto vals = harmonicFaces[j].evaluateAtPoints(globalToLocalIdMap(i, j), faceSamples2d[j]);
            for(auto x : vals) b(cnt++) = x;
        }
        
        // now evaluate at edge/vertex probes
        b.bottomRows(PEdgeSamples.rows()) = PEdgeSamples * Eigen::VectorXd::Unit((int)vertices.size(), i);
        
        // solve for coefficients (one for each rbf, 4 for the linear part)
        coefficients.row(i) = chol.solve(A.transpose() * b);
    }
}

