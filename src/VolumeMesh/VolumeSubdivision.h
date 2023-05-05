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
    explicit VolumeSubdivision(VolumeMesh& mesh);

    void tetrahedra(int face_point, int cell_point);

    void irregular_mesh(int n);

    void full_truncation();

    void quads();

    void linear_subdivision();

private:
    VolumeMesh& mesh_;
};
