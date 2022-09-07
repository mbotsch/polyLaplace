//=============================================================================

#include "SurfaceViewer.h"

//=============================================================================
int main(int argc, char **argv) {
    if(argc < 2){
        std::cout << "Specify input mesh" << std::endl;
    }
    Viewer window("Polygon Modeling", 800, 600);
    window.load_mesh(argv[1]);
    return window.run();
}
//=============================================================================
