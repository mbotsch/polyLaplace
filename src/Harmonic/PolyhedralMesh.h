#pragma once

#include <vector>
#include <array>
#include <string>
#include <Eigen/Dense>

class PolyhedralMesh
{
public:
    typedef int Index;
    typedef std::pair<Index, char> HalfIndex;
    typedef std::array<double, 3> Point;
    typedef std::array<Index, 2> Edge;

    std::vector<Point> points;
    std::vector<Edge> edges;
    std::vector<std::vector<HalfIndex>> faces; // list of half edges
    std::vector<std::vector<HalfIndex>> cells; // list of half faces
    std::vector<char> isBoundary;

    // return list of face vertex indizes
    std::vector<Index> faceIndizes(const HalfIndex& halfFace) const;

    // return list of vertex index lists
    std::vector<std::vector<Index>> cellfaceIndizes(
        const Index cellIndex) const;

    bool loadOVM(const std::string& filename);

    void identifyBoundary();

public:
    PolyhedralMesh();

    PolyhedralMesh(const std::string& filename);

    // get surface mesh for a single cell
    std::vector<Index> getCellGeometry(const Index cellIndex,
                                       std::vector<Point>& pts,
                                       std::vector<std::vector<Index>>& fcs,
                                       bool triangularize = false) const;

    std::vector<Index> getCellGeometry(
        const Index cellIndex, std::vector<Eigen::Vector3d>& cellPoints,
        std::vector<std::vector<Index>>& cellFaces,
        bool triangularize = false) const;

    int numCells() const;

    int numVertices() const;

    bool isBoundaryVertex(const int i) const;
};
