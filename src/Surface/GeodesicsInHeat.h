
//=============================================================================
#pragma once
//=============================================================================
#define PMP_SCALAR_TYPE_64

 #include <pmp/SurfaceMesh.h>
#include "PolyLaplace.h"
#include <Eigen/Sparse>
#include <fstream>
#include <pmp/algorithms/SurfaceNormals.h>
//=============================================================================

class GeodesicsInHeat
{
public:
    GeodesicsInHeat(pmp::SurfaceMesh& mesh, int laplace, int min_point,
                    bool geodist, bool euklid, bool mean_edge_ = false);
    ~GeodesicsInHeat();

    void getDistance(const int vertex, Eigen::VectorXd& dist,
                     Eigen::VectorXd& orthodist);

    void distance_to_texture_coordinates() const;

    void compute_geodesics(bool lumped = true);

    void write_Distance_error(const int vertex, Eigen::VectorXd& dist,
                     Eigen::VectorXd& orthodist, bool lumped,std::ofstream& file);

private:
    pmp::SurfaceMesh& mesh_;

    unsigned int laplace_, min_point_;

    Eigen::MatrixXd pos;

    Eigen::SparseMatrix<double> divOperator, gradOperator;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cholL, cholA;

    double averageEdgeLength(const pmp::SurfaceMesh& mesh);

    double maxEdgeLength(const pmp::SurfaceMesh& mesh);

    void buildDivOperator();

    void buildGradientOperator();

    double great_circle_distance(pmp::Vertex v, pmp::Vertex vv, double r = 1.0);

    double haversine_distance(pmp::Vertex v, pmp::Vertex vv, double r = 1.0);

    double vincenty_distance(pmp::Vertex v, pmp::Vertex vv, double r = 1.0);

    bool geodist_sphere_, geodist_cube_;

    bool mean_edge_ ;

};
