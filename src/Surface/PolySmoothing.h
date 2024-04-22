//=============================================================================
// Copyright 2024 Astrid Bunge, Sven Wagner, Dennis Bukenberger, Mario Botsch, Marc Alexa
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include <pmp/surface_mesh.h>
#include "diffgeo.h"
#include "LaplaceConstruction.h"
#include "pmp/bounding_box.h"

//=============================================================================

using namespace pmp;

class PolySmoothing
{
public:
    PolySmoothing(SurfaceMesh& mesh, SmoothingConfigs oConf);
    void optimize(int quadricsTau = 5);
    void computeVirtualVertices(bool use_fallback = true);
    void setupQuadrics();

private:
    SurfaceMesh& mesh_;
    SmoothingConfigs oConf_;
    bool is2D_;
    int penalty_;
};