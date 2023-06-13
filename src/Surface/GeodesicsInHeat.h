//=============================================================================
// Copyright 2023 Astrid Bunge, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================
#include <pmp/SurfaceMesh.h>
#include "[AW11]Laplace.h"
#include <Eigen/Sparse>
//=============================================================================
enum DiffusionStep
{
    MeanEdge = 0,
    MaxEdge = 1,
    MaxDiagonal = 2
};
class GeodesicsInHeat
{
public:
    GeodesicsInHeat(pmp::SurfaceMesh& mesh, int laplace, int min_point,
                    bool geodist, bool euklid,
                    DiffusionStep diffusion = MeanEdge);
    ~GeodesicsInHeat();

    double getDistance(int vertex, Eigen::VectorXd& dist,
                       Eigen::VectorXd& orthodist);

    void distance_to_texture_coordinates() const;

    void compute_geodesics(bool lumped = true);

private:
    pmp::SurfaceMesh& mesh_;

    int laplace_, min_point_;

    Eigen::MatrixXd pos;

    Eigen::SparseMatrix<double> divOperator, gradOperator;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cholL, cholA;

    static double edgeLength(const pmp::SurfaceMesh& mesh, pmp::Edge e);
    static double averageEdgeLength(const pmp::SurfaceMesh& mesh);
    static double maxEdgeLength(const pmp::SurfaceMesh& mesh);
    static double maxDiagonalLength(const pmp::SurfaceMesh& mesh);

    double great_circle_distance(pmp::Vertex v, pmp::Vertex vv, double r = 1.0);
    double haversine_distance(pmp::Vertex v, pmp::Vertex vv, double r = 1.0);
    double vincenty_distance(pmp::Vertex v, pmp::Vertex vv, double r = 1.0);

    bool geodist_sphere_, geodist_cube_;

    DiffusionStep diffusionStep_;
};
