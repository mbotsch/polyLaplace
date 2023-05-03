# PolygonLaplacians

Since we use [pmp-library](http://www.pmp-library.org/) as submodule, you have to clone the repository recursively:

    git clone --recursive git@github.com:mbotsch/polyLaplace.git

Configure and build:

    cd polyLaplace && mkdir build && cd build && cmake .. && make

This will automatically build our code and all dependencies. Finally, start the GUI app with a polygon mesh or a polyhedral mesh:

    ./surface_viewer ../data/surface_meshes/grid/quad_4.obj 
    ./volume_viewer ../data/volume_meshes/cubes/cube_hexahedron_3.ovm 

Alternatively you can either run the volume and surface convergence tests by executing:

    ./run_surface_tests
    ./run_volume_tests

Have fun!
