//=============================================================================

#include "VolumeViewer.h"

//=============================================================================


int main(int argc, char **argv) {

    VolumeViewer volumeWindow("Polygon Modeling", 800, 600);
    if(argc < 2){
        volumeWindow.load_mesh("../data/ovm/cube.ovm");
    }else{
//        for(int i = 1; i< argc; i++){
//            volumeWindow.write_histogram(argv[i]);
//        }
        volumeWindow.load_mesh(argv[1]);
    }
    return volumeWindow.run();

}


//=============================================================================
