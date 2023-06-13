//=============================================================================
// Copyright 2023 Astrid Bunge, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include <pmp/SurfaceMesh.h>

//=============================================================================

using namespace pmp;

//=============================================================================

class Parameterization
{
public:
    // give a mesh in the constructor
    Parameterization(SurfaceMesh& mesh) : mesh_(mesh) {}

    // compute discrete harmonic parameterization
    bool harmonic(unsigned int laplace, unsigned int min_point);

private:
    // map boundary to unit circle
    bool setup_boundary_constraints();

private:
    SurfaceMesh& mesh_;
};

//=============================================================================
