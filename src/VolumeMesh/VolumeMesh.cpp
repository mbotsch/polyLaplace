//=============================================================================
// Copyright 2021 Hendrik Meyer, Astrid Bunge, Mario Botsch
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "VolumeMesh.h"
#include "VolumeMeshIO.h"

//=============================================================================

VolumeMesh::VolumeMesh() = default;

//-----------------------------------------------------------------------------

VolumeMesh::~VolumeMesh() = default;

//-----------------------------------------------------------------------------

bool VolumeMesh::read(const std::string &filename)
{
    VolumeMeshIO reader(filename);
    return reader.read(*this);
}

//-----------------------------------------------------------------------------

void VolumeMesh::update_bounding_sphere()
{
    Vec3d v;
    Vec3d min(std::numeric_limits<double>::max());
    Vec3d max(std::numeric_limits<double>::min());

    for (auto bv_it = bv_iter(); bv_it.valid(); ++bv_it)
    {
        v = vertex(*bv_it);
        min.minimize(v);
        max.maximize(v);
    }

    Vec3d center, diameter;
    double radius;

    center = 0.5 * (min + max);
    diameter = max - min;
    radius = 0.5 * diameter.norm();

    bounding_sphere_ = BoundingSphere(center, radius);
}
