# Discrete Laplacians for  Polygonal/Polyhedral Meshes

This repository contains the source code for the discrete polygonal/polyhedral Laplacians proposed in these papers:

- Bunge, Herholz, Kazhdan, Botsch,
  **Polygon Laplacian Made Simple**,
  Computer Graphics Forum 39(2), 2020.
- Bunge, Botsch, Alexa,
  **The Diamond Laplace for Polygonal and Polyhedral Meshes**,
  Computer Graphics Forum 40(5), 2021.
- Bunge, Bukenberger, Wagner, Alexa, Botsch,
  **Polygon Laplacian Made Robust**,
  Computer Graphics Forum 43(2), 2024.  

In order to compare to related methods, the code also contains implementations of these methods:

- Martin, Kaufmann, Botsch, Wicke, Gross,
  **Polyhedral Finite Elements Using Harmonic Basis Functions**,
  Computer Graphics Forum 27(5), 2008.
- Alexa, Wardetzky,
  **Discrete Laplacians on General Polygonal Meshes**,
  ACM Transactions on Graphics 30(4), 2011.
- de Goes, Butts, Desbrun,
  **Discrete Differential Operators on Polygonal Meshes**,
  ACM Transactions on Graphics 39(4), 2020.

Comparisons of the different Laplace operators on a range of example applications (also included in the source code) are provided in these papers:

- Bunge, Botsch,
  **A Survey on Discrete Laplacians for General Polygonal Meshes**,
  Computer Graphics Forum 42(2), 2023.
- Bunge, Alexa, Botsch,
  **Discrete Laplacians for General Polygonal and Polyhedral Meshes**,
  SIGGRAPH Asia Courses, 2023.


## Installation, building, and running the demos

Since we use [pmp-library](http://www.pmp-library.org/) as submodule, you have to clone the repository recursively:

    git clone --recursive git@github.com:mbotsch/polyLaplace.git

Configure and build:

    cd polyLaplace && mkdir build && cd build && cmake .. && make

This will automatically build our code and all dependencies. Finally, start the GUI apps with a polygon mesh or a polyhedral mesh:

    ./poly_demo ../data/fertility.obj
    ./surface_viewer ../data/surface_meshes/grid/quad_4.obj 
    ./volume_viewer ../data/volume_meshes/cubes/cube_hexahedron_3.ovm 

Alternatively you can run the volume and surface convergence tests by executing:

    ./run_surface_tests
    ./run_volume_tests

Have fun!
