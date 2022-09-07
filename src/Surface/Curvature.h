//=============================================================================
#pragma once
//=============================================================================

#include <pmp/SurfaceMesh.h>
#include "diffgeo.h"
#include "[AW11]Laplace.h"
#include "LaplaceConstruction.h"
#include "[dGBD20]Laplace.h"
#include <pmp/algorithms/Normals.h>
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
    void visualize_curvature(int laplace, int min_point_, bool lumped=true);

    double compute_curvature_error(int laplace,  int min_point_, bool lumped=true);

private:
    SurfaceMesh& mesh_;
    bool compare_to_sphere;

    //! convert curvature values ("v:curv") to 1D texture coordinates
    void curvature_to_texture_coordinates() const;
};

//=============================================================================
