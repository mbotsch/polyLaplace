# Installation {#installation}

In this section, we describe how to configure, build, and install the
pmp-library in detail.

## System Requirements

The pmp-library uses [CMake](http://www.cmake.org) as its build configuration
system. Version 3.16.3 or greater is required. The pmp-library requires a
C++14-compliant compiler. We continuously build and test the pmp-library
with the following compilers and operating systems:

| Operating System | Compiler           |
| ---------------- | ------------------ |
| Linux            | gcc 9.3.0          |
| macOS            | AppleClang 12.0.0  |
| Windows          | Visual Studio 2019 |

## Dependencies

Some parts of the pmp-library depend on the following third-party libraries:

| Library                                             | Description                       | Version     |
| --------------------------------------------------- | --------------------------------- | ----------- |
| [Eigen](http://eigen.tuxfamily.org)                 | C++ linear algebra library        | &ge; 3.4.0  |
| [OpenGL](http://opengl.org)                         | Open Graphics Library             | &ge; 3.3    |
| [GLEW](http://glew.sourceforge.net)                 | OpenGL Extension Wrangler Library | &ge; 2.1.0  |
| [GLFW](http://glfw.org)                             | Graphics Library Framework        | &ge; 3.3.8  |
| [ImGui](https://github.com/ocornut/imgui)           | Immediate Mode GUI                | &ge; 1.70   |
| [Google Test](https://github.com/google/googletest) | C++ Test Framework                | &ge; 1.12.1 |

By default, we include the corresponding libraries using git submodules. Note
that OpenGL and related dependencies are optional. They are only needed if you
want to use the viewer classes. Google Test is optional as well and only
required if you want to run the unit test suite.

## Configuration

The pmp-library relies on [CMake](http://www.cmake.org) as its build and
configuration system. `CMake` is a cross-platform build-system capable of
generating different build files (so-called _generators_) for a specific
platform, e.g., Makefiles for Linux/Unix, Xcode projects for Mac OS-X and Visual
Studio projects for Windows.

On the command line change to the top-level pmp-library directory, create a
build directory and run `cmake`:

```sh
cd pmp-library
mkdir build
cd build
cmake ..
```

The configuration procedure can be fine-tuned by specifying flags using the `-D`
option of `cmake`:

```sh
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/g++
```

The command above would configure `CMake` to use release mode as its build type
and `/usr/bin/g++` as its C++ compiler.

In order to compile the included examples configure `CMake` with

```sh
cmake -DWITH_EXAMPLES=true ..
```

Commonly used flags are shown below.

| Flag                 | Description                                        |
| -------------------- | -------------------------------------------------- |
| `CMAKE_BUILD_TYPE`   | Specify the build type, e.g. Debug or Release.     |
| `CMAKE_CXX_COMPILER` | Specify the compiler to be used.                   |
| `CMAKE_CXX_FLAGS`    | Specify additional compiler flags, e.g. `-DNDEBUG` |

For additional information on using `CMake` and
customizing its configuration see
the [CMake documentation](http://cmake.org/cmake/help/documentation.html).

## Building

After successful configuration pmp-library can be build using the chosen build
system. For a Unix-like environment the default generator is Makefiles. In order
to build pmp-library just call

```sh
make
```

from the top-level build directory. In order to build pmp in parallel use the
`-j` option of `make`:

```sh
make -j
```

The resulting library is named `libpmp.so` and
located in the current working directory.

In order to build the full HTML manual and reference documentation call

```sh
make docs
```

The resulting HTML documentation can be found in the `docs/html/` sub-directory.
Note: this requires [Doxygen](http://www.doxygen.nl/) to be installed. In order
to generate proper bibliographical references please install
[BibTex](http://www.bibtex.org/) as well.

## Installation

In order to install pmp-library just call

```sh
sudo make install
```

Upon installation, both the library and headers will be installed to the
directory given via `CMAKE_INSTALL_PREFIX`, which defaults to `/usr/local/` on
Unix-like systems. If you need to install to a custom location set the install
prefix during build configuration:

```sh
cmake -DCMAKE_INSTALL_PREFIX=<your custom path> ..
```

The library can be uninstalled using

```sh
make uninstall
```

To use the pmp-library in your own CMake-based projects simply include the
library by using `find_package(pmp)` and point CMake to the directory containing
the pmp-library CMake configuration file `pmpConfig.cmake`. This can be either
the pmp-library build directory

```sh
cmake -Dpmp_DIR=<path-to-pmp-build-directory>
```

or the installed version

```sh
cmake -Dpmp_DIR=<your custom path>/lib/cmake/pmp
```

This way, you can simply link your own target against pmp-library

```cmake
target_link_libraries(your_target pmp)
```

**Note:** The usage described above is currently limited to the @ref core and
@ref algorithms modules of the pmp-library. If you want to use the @ref
visualization module you need to link your target against `pmp_vis` and its
dependencies: `stb_image`, `imgui`, `glfw`, `glew`, as well as your platform
OpenGL library.

## Build Options

### Index Type

By default, the pmp-libray uses 32-bit unsigned integers as internal index type
to reference entities. However, if you need to process very large data sets this
might not be sufficient. In this case, you can change the index type to be
64-bit by specifying

```sh
cmake -DPMP_INDEX_TYPE=64
```

during build configuration.

### Scalar Type

By default, the pmp-library uses `float` as `Scalar` type. In case you require
higher floating point precision you can change the `Scalar` type to `double` by
specifying

```sh
cmake -DPMP_SCALAR_TYPE=64
```

during build configuration.

## Building JavaScript Apps

In order to build the JavaScript/WebAssembly applications
using [emscripten](https://github.com/kripken/emscripten), download the SDK
from <https://github.com/kripken/emscripten> and follow the installation
instructions.

Next, source the environment setup script:

```sh
source <path_to_install_dir>/emsdk_env.sh
```

Create a build directory, run cmake and build:

```sh
mkdir html
cd html
emcmake cmake ..
make
```

Finally, start a local webserver and open the HTML apps:

```sh
python3 -m http.server
<your-browser> localhost:8000
```

You can also run HTML/WASM apps using `emrun mpview.html`, but you might have to adjust
some [browser settings](https://emscripten.org/docs/compiling/Running-html-files-with-emrun.html) first.
