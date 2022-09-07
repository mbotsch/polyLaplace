// Copyright 2021 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#pragma once

#include <pmp/SurfaceMesh.h>

namespace pmp {

// generate triangle fan of six triangles around center vertex
SurfaceMesh vertex_onering();

// generate onering around an edge
SurfaceMesh edge_onering();

// generate subdivided icosahedron
// based on Loop subdivision
// original icosahedron edges are marked as feature edges
SurfaceMesh subdivided_icosahedron();

// generate 2d non-convex L shape
SurfaceMesh l_shape();

// generate cone with bottom polygon removed
SurfaceMesh open_cone();

// a mesh with texcoords and texture seams
SurfaceMesh texture_seams_mesh();

} // namespace pmp