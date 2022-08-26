//=============================================================================
#pragma once
//=============================================================================

#include <pmp/SurfaceMesh.h>
#include "diffgeo.h"
#include "PolyLaplace.h"
#include "AQAPoly_Laplacian.h"
#define PMP_SCALAR_TYPE_64

//=============================================================================

using namespace pmp;

//=============================================================================

class Curvature
{
public:
    Curvature(SurfaceMesh& mesh, bool compare)
            : mesh_(mesh), compare_to_sphere(compare)
    {
    }

    //! Visualizes the mean curvature of our mesh.
    void visualize_curvature(unsigned int laplace, unsigned int min_point_, bool lumped=true, int degree = 1,CoarseDimension coarseningType=Edges);

    double compute_curvature_error(unsigned int laplace, unsigned int min_point_, bool lumped=true, int degree = 1,CoarseDimension coarseningType=Edges);

    double compute_quad_Tri_Laplace_curvature_error(SurfaceMesh& initial_mesh);
private:
    SurfaceMesh& mesh_;
    bool compare_to_sphere;

    //! convert curvature values ("v:curv") to 1D texture coordinates
    void curvature_to_texture_coordinates() const;
};

//=============================================================================
