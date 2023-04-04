//=============================================================================

#include "SurfaceViewer.h"

//=============================================================================
int main(int argc, char **argv) {
    Viewer window("Polygon Modeling", 800, 600);
    if(argc < 2){
        window.load_mesh("../data/surface_meshes/grid/concave_1.obj");
    }else {
        window.load_mesh(argv[1]);
    }
    return window.run();
}
//=============================================================================
