//=============================================================================
#pragma once
//=============================================================================

#include <pmp/SurfaceMesh.h>
#include "diffgeo.h"
#include "LaplaceConstruction.h"
#include "pmp/BoundingBox.h"

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