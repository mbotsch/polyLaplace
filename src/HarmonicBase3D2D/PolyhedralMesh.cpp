#include "PolyhedralMesh.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>
#include <map>

using namespace std;

PolyhedralMesh::PolyhedralMesh() {}

PolyhedralMesh::PolyhedralMesh(const std::string filename) {
    loadOVM(filename);
    identifyBoundary();
}

int PolyhedralMesh::numCells() const {
    return (int)cells.size();
}

int PolyhedralMesh::numVertices() const {
    return (int)points.size();
}


bool PolyhedralMesh::isBoundaryVertex(const int i) const {
    return isBoundary[i];
}

void PolyhedralMesh::identifyBoundary() {
    std::map<int, char> faceCount;
    
    for(auto& c : cells) {
        for(auto& f : c) {
            ++faceCount[f.first];
        }
    }
    
    isBoundary.clear();
    isBoundary.resize(points.size());
    
    for(auto& fc : faceCount) {
        if(fc.second == 1) {
            for(auto f : faces[fc.first]) {
                isBoundary[edges[f.first][0]] = 1;
                isBoundary[edges[f.first][1]] = 1;
            }
        }
    }
}

bool PolyhedralMesh::loadOVM(const std::string filename) {
    
    ifstream file(filename);
    if(!file.is_open()) return false;
    
    auto CheckString = [&](const std::string& s) {
        std::string buf;
        file >> buf;
        if(s.compare(buf)) {
            std::cout << "reading ovm: string " << s << " not found. Returning." << std::endl;
            return false;
        }
        
        return true;
    };
    
    if(!CheckString("OVM") || !CheckString("ASCII") || !CheckString("Vertices")) return false;
    
    int nv;
    file >> nv;
    points.reserve(nv);
    
    double x, y, z;
    for(int i = 0; i < nv; ++i) {
        file >> x >> y >> z;
        points.push_back({x, y, z});
    }
    
    if(!CheckString("Edges")) return false;
    int ne;
    file >> ne;
    edges.reserve(ne);
    
    Edge e;
    for(int i = 0; i < ne; ++i) {
        file >> e[0] >> e[1];
        edges.push_back(e);
    }
    
    if(!CheckString("Faces")) return false;
    int nf;
    file >> nf;
    
    for(int i = 0; i < nf; ++i) {
        int nei;
        file >> nei;
        faces.reserve(nei);
        
        std::vector<HalfIndex> face(nei);
        Index id;
        for(int j = 0; j < nei; ++j) {
            file >> id;
            face[j] = make_pair(id >> 1, (char)(id & 1));
        }
        
        faces.push_back(move(face));
    }
    
    if(!CheckString("Polyhedra")) return false;
    int nc;
    file >> nc;
    cells.reserve(nc);
    
    for(int i = 0; i < nc; ++i) {
        int nfi;
        file >> nfi;
        std::vector<HalfIndex> cell(nfi);
        Index id;
        for(int j = 0; j < nfi; ++j) {
            file >> id;
            cell[j] = make_pair(id >> 1, (char)(id & 1));
        }
        
        cells.push_back(move(cell));
    }
    
    file.close();
    
    return true;
}

std::vector<PolyhedralMesh::Index> PolyhedralMesh::getCellGeometry(const Index cellIndex, std::vector<Eigen::Vector3d>& cellPoints, std::vector<std::vector<Index>>& cellFaces, bool triangularize) const {
    
    std::vector<Point> pts;
    auto ret = getCellGeometry(cellIndex, pts, cellFaces, triangularize);
    
    for(auto& p : pts) cellPoints.push_back(Eigen::Vector3d(p[0], p[1], p[2]));
    
    return ret;
}

std::vector<PolyhedralMesh::Index> PolyhedralMesh::getCellGeometry(const Index cellIndex, std::vector<Point>& cellPoints, std::vector<std::vector<Index>>& cellFaces, bool triangularize) const {
    
    auto faces = cellfaceIndizes(cellIndex);
    
    vector<Index> vertices;
    for(auto& f : faces) copy(f.begin(), f.end(), back_inserter(vertices));
    sort(vertices.begin(), vertices.end());
    vertices.erase(unique(vertices.begin(), vertices.end()), vertices.end());
    
    map<Index, Index> vertexMap;
    for(int i = 0; i < vertices.size(); ++i) vertexMap[vertices[i]] = i;
    
    for(auto i : vertices) cellPoints.push_back(points[i]);
    
    
    for(auto& f : faces) {
        vector<Index> polyVerts;
        
        for(int k : f) polyVerts.push_back(vertexMap[k]);
        reverse(polyVerts.begin(), polyVerts.end());
        
        if(triangularize) {
            for(int i = 1; i < polyVerts.size() - 1; ++i) {
                cellFaces.push_back({polyVerts[0],  polyVerts[i],  polyVerts[i+1]});
            }
        } else cellFaces.push_back(polyVerts);
    }
    
    return vertices;
}
  
vector<PolyhedralMesh::Index> PolyhedralMesh::faceIndizes(const HalfIndex& halfFace) const {
    std::vector<Index> ids;
    
    for(const auto& j : faces[halfFace.first]) {
        ids.push_back(j.second ? edges[j.first][1] : edges[j.first][0]);
    }
    
    if(halfFace.second) reverse(ids.begin(), ids.end());
    
    return ids;
}

vector<vector<PolyhedralMesh::Index>> PolyhedralMesh::cellfaceIndizes(const Index cellIndex) const {
    vector<vector<Index>> ids;
    
    for(const auto& hf : cells[cellIndex]) {
        ids.push_back(faceIndizes(hf));
    }

    return ids;
}

