
//=============================================================================
#pragma once
//=============================================================================
#define PMP_SCALAR_TYPE_64

 #include <pmp/SurfaceMesh.h>
#include "[AW11]Laplace.h"
#include <Eigen/Sparse>
#include <pmp/algorithms/Normals.h>
//=============================================================================

class GeodesicsInHeat
{
public:
    GeodesicsInHeat(pmp::SurfaceMesh& mesh, int laplace, int min_point,
                    bool geodist, bool euklid, bool mean_edge_ = false);
    ~GeodesicsInHeat();

    void getDistance(int vertex, Eigen::VectorXd& dist,
                     Eigen::VectorXd& orthodist);

    void distance_to_texture_coordinates() const;

    void compute_geodesics(bool lumped = true);

private:
    pmp::SurfaceMesh& mesh_;

    int laplace_, min_point_;

    Eigen::MatrixXd pos;

    Eigen::SparseMatrix<double> divOperator, gradOperator;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cholL, cholA;

    static double averageEdgeLength(const pmp::SurfaceMesh& mesh);

    static double maxEdgeLength(const pmp::SurfaceMesh& mesh);

    double great_circle_distance(pmp::Vertex v, pmp::Vertex vv, double r = 1.0);

    double haversine_distance(pmp::Vertex v, pmp::Vertex vv, double r = 1.0);

    double vincenty_distance(pmp::Vertex v, pmp::Vertex vv, double r = 1.0);

    bool geodist_sphere_, geodist_cube_;

    bool mean_edge_ ;

};
