#include "WriteTetOVM.hpp"
// Include the file manager header
#include <OpenVolumeMesh/FileManager/FileManager.hh>
// Include the polyhedral mesh header
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>

#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>

void writeTetOVM(const std::vector<Eigen::Vector3d>& pts, const std::vector<std::vector<int>>& tets, const std::string& fname) {

    using namespace OpenVolumeMesh;
    
    GeometricTetrahedralMeshV3d mesh;
    
    std::map<int, VertexHandle> vmap;
    
    // add vertices
    for(int i = 0; i < (int)pts.size(); ++i) {
        auto p = pts[i];
        vmap[i] = mesh.add_vertex(Vec3d(p(0), p(1), p(2)));
    }
    
    // add cells
    for(auto& t : tets) {
        mesh.add_cell(vmap[t[0]], vmap[t[1]], vmap[t[2]], vmap[t[3]], true);
    };
    
    OpenVolumeMesh::IO::FileManager fileManager;
    fileManager.writeFile(fname, mesh);
}
