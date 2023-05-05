//=============================================================================

#include "VolumeViewer.h"

//=============================================================================

int main(int argc, char** argv)
{
    VolumeViewer volumeWindow("Polygon Modeling", 800, 600);
    if (argc < 2)
    {
        volumeWindow.load_mesh(
            "../data/volume_meshes/cubes/cube_hexahedron_1.ovm");
    }
    else
    {
        volumeWindow.load_mesh(argv[1]);
    }
    return volumeWindow.run();
}

//=============================================================================
