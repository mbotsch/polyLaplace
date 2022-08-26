//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "SurfaceViewer.h"
#include "VolumeViewer.h"

//=============================================================================

int choose_viewer(std::string filename)
{
    std::string::size_type dot(filename.rfind("."));
    if (dot == std::string::npos)
        return false;
    std::string ext = filename.substr(dot + 1, filename.length() - dot - 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), tolower);

    int viewer = 0;

    if (ext == "ovm" || ext == "mesh")
    {
        viewer = 1;
    }

    return viewer;
}

int main(int argc, char **argv)
{
    int viewer = 0;
    if (argc == 2)
    {
        viewer = choose_viewer(argv[1]);
    }

    if (viewer == 0)
    {
        Viewer window("Polygon Modeling", 800, 600);
        window.load_mesh(argv[1]);
        return window.run();
    }
    else
    {
        VolumeViewer volumeWindow("Polygon Modeling", 800, 600);
        volumeWindow.load_mesh(argv[1]);
        return volumeWindow.run();
    }
}

//=============================================================================
