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
//typedef OpenVolumeMesh::Vec3f Vec3f;
typedef OpenVolumeMesh::Vec3d Vec3d;

//typedef OpenVolumeMesh::GeometricPolyhedralMeshV3f PolyhedralMeshV3f;
typedef OpenVolumeMesh::GeometricPolyhedralMeshV3d PolyhedralMeshV3d;

typedef OpenVolumeMesh::HalfFaceHandle HFHandle;
typedef OpenVolumeMesh::VertexHandle VHandle;
typedef OpenVolumeMesh::CellHandle CHandle;

//typedef Eigen::Matrix<float, Eigen::Dynamic, 3> VertVec;
//typedef Eigen::Triplet<float> T;
//typedef Eigen::SparseMatrix<float> SpMat;

typedef Eigen::Matrix<double, Eigen::Dynamic, 3> VertVec;
typedef Eigen::Triplet<double> T;
typedef Eigen::SparseMatrix<double> SpMat;

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

    static bool write_cell_histogram(const std::string &filename);

public:
    struct BoundingSphere
    {
        //        explicit BoundingSphere(Vec3f center = Vec3f(0, 0, 0), float radius = 1)
        //            : center(center), radius(radius)
        //        {
        //        }
        //
        //        Vec3f center;
        //        float radius;

        explicit BoundingSphere(Vec3d center = Vec3d(0, 0, 0),
                                double radius = 1)
            : center(center), radius(radius)
        {
        }

        Vec3d center;
        double radius;
    };

    //! return bounding sphere enclosig the mesh
    //    std::pair<Vec3f, float> get_bounding_sphere()
    //    {
    //        return {bounding_sphere_.center, bounding_sphere_.radius};
    //    }
    std::pair<Vec3d, double> get_bounding_sphere()
    {
        return {bounding_sphere_.center, bounding_sphere_.radius};
    }

    //! calculate bounding sphere enclosing the mesh
    void update_bounding_sphere();

private:
    BoundingSphere bounding_sphere_; //! bounding sphere enclosing the mesh
};
