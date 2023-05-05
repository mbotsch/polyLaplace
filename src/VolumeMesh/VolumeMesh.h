//=============================================================================
// Copyright 2021 Hendrik Meyer, Astrid Bunge, Mario Botsch
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <Eigen/SparseCore>

//=============================================================================

#define OVM OpenVolumeMesh

typedef OpenVolumeMesh::Vec3d Vec3d;

typedef OpenVolumeMesh::GeometricPolyhedralMeshV3d PolyhedralMeshV3d;

typedef OpenVolumeMesh::HalfFaceHandle HFHandle;
typedef OpenVolumeMesh::VertexHandle VHandle;
typedef OpenVolumeMesh::CellHandle CHandle;

typedef Eigen::Triplet<double> T;

//=============================================================================

class VolumeMesh : public PolyhedralMeshV3d
{
public:
    //! constructor
    VolumeMesh();

    //! destructor
    ~VolumeMesh() override;

    //! read file
    bool read(const std::string& filename);

public:
    struct BoundingSphere
    {
        explicit BoundingSphere(Vec3d center = Vec3d(0, 0, 0),
                                double radius = 1)
            : center(center), radius(radius)
        {
        }

        Vec3d center;
        double radius;
    };

    //! return bounding sphere enclosig the mesh
    std::pair<Vec3d, double> get_bounding_sphere()
    {
        return {bounding_sphere_.center, bounding_sphere_.radius};
    }

    //! calculate bounding sphere enclosing the mesh
    void update_bounding_sphere();

private:
    BoundingSphere bounding_sphere_; //! bounding sphere enclosing the mesh
};
