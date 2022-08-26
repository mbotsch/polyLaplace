//=============================================================================
// Copyright 2021 Hendrik Meyer, Astrid Bunge, Mario Botsch
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include "VolumeMesh.h"

//=============================================================================

class VolumeSubdivision
{
public:
    VolumeSubdivision(VolumeMesh& mesh);

    void tetrahedra(unsigned int face_point, unsigned int cell_point);

    void irregular_mesh(int n);

    void full_truncation();

    void quads();

    void linear_subdivision();

private:
    VolumeMesh& mesh_;
};
