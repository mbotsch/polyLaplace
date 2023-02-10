#pragma once

#include "GaussQuadrature.h"
#include "HarmonicPolygon.hpp"
#include "Geometry.hpp"
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fstream>

struct HEdge {
    
    int i;
    int j;
    
    HEdge(const int i_, const int j_);
    
    bool operator==(const HEdge& e2) const;
    
    bool operator<(const HEdge& e2) const;
};

struct Normal {
    Eigen::Vector3d n;
    
    Normal();
    
    void normalize();
    
    const Normal& operator+=(const Eigen::Vector3d& v);
};
   

void loadPoly(const std::string& fname, std::vector<Eigen::Vector3d>& vertices, std::vector<std::vector<int>>& faces);

class HarmonicPolyhedron {
public:

    double eps = 1e-2;

    int probesPerEdge = 9;
    int probesPerFace = 30;
    int kernelsPerFace = 10;
    int kernelsPerEdge = 3;

    double scale{};
    Eigen::Vector3d mean;
    
    TetrahedralQuadrature tetQuad;
    
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::vector<int>> faces;
    std::vector<std::vector<int>> tris;
    std::vector<std::vector<Eigen::Vector3d>> triPoints;
    
    Eigen::MatrixXd coefficients;
    
    int nProbes = 0;
    int nKernels = 0;
    
    std::vector<HarmonicPolygon> harmonicFaces;
    std::vector<std::vector<Eigen::Vector2d>> faceSamples2d;
    std::vector<Eigen::Vector3d> kernelCenters;
    
    std::vector<HEdge> edges;
    Eigen::SparseMatrix<double> PEdgeSamples; // linearly interpolates values per vertex to edge samples
    Eigen::MatrixXd edgeSamples;
    
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
    void dumpProbesAndKernels();
    
    void dumpCrossection();
    
    bool insideConvex(const Eigen::Vector3d& p);
    
    bool inside(const Eigen::Vector3d& p);
    
    int globalToLocalIdMap(const int vid, const int fid);
    
    Eigen::Vector3d evaluateGrad(const int node, const Eigen::Vector3d& p);
    
    double evaluate(const int node, const Eigen::Vector3d& p);
    
    void initQuadrature();
    
    void stiffnessMatrix(Eigen::MatrixXd& K);
    
    void massMatrix(Eigen::MatrixXd& M);
    
    HarmonicPolyhedron();
    
    HarmonicPolyhedron(const std::string& fname);
    
    HarmonicPolyhedron(std::vector<Eigen::Vector3d>  vertices_, std::vector<std::vector<int>>  faces_);
    
    void init();
};
