//=============================================================================
// Copyright 2021 Hendrik Meyer, Astrid Bunge, Mario Botsch
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include "VolumeMesh.h"

//=============================================================================

class VolumeMeshIO
{
public:
    VolumeMeshIO(const std::string& filename) : filename_(filename) {}

    bool read(VolumeMesh& mesh);

    bool write_histrogram();

    bool write(const VolumeMesh& mesh);

private:
    bool read_ovm(VolumeMesh& mesh);
    bool read_mesh(VolumeMesh& mesh);
    bool read_hybrid(VolumeMesh& mesh);

    bool write_ovm(const VolumeMesh& mesh);
    //    bool write_mesh(const VolumeMesh& mesh);
    bool write_hybrid(const VolumeMesh& mesh);

private:
    std::string filename_;
};
