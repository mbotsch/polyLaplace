project(PolygonLaplacians)

cmake_minimum_required(VERSION 3.0.2)

# project version
set(CMAKE_PROJECT_VERSION 1.0.0)
set(CMAKE_PROJECT_VERSION_MAJOR 1)
set(CMAKE_PROJECT_VERSION_MINOR 0)
set(CMAKE_PROJECT_VERSION_PATCH 0)

# set output paths
set(PROJECT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

# set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# compile PMP library
set(PMP_BUILD_APPS     OFF CACHE BOOL "")
set(PMP_BUILD_EXAMPLES OFF CACHE BOOL "")
set(PMP_BUILD_TESTS    OFF CACHE BOOL "")
set(PMP_BUILD_DOCS     OFF CACHE BOOL "")
#set(PMP_SCALAR_TYPE 64)

add_subdirectory(external/pmp-library)

add_subdirectory(external/OpenVolumeMesh)

# add include directories
include_directories(${PROJECT_SOURCE_DIR}/external/pmp-library/src)
include_directories(${PROJECT_SOURCE_DIR}/external/pmp-library/external/imgui)
include_directories(${PROJECT_SOURCE_DIR}/external/pmp-library/external/eigen)
include_directories(${PROJECT_SOURCE_DIR}/external/libigl/include)
include_directories(${PROJECT_SOURCE_DIR}/external/spectra/include)
include_directories(${PROJECT_SOURCE_DIR}/external/OpenVolumeMesh/src)

#Path to cholmod
#include_directories(/home/abunge/work/SuiteSparse-5.8.1/CHOLMOD/Include/)

#Path to triangle
#include_directories(/home/abunge/programs/libigl/external/triangle/)

add_definitions(-DINCLUDE_TEMPLATES)

if (NOT EMSCRIPTEN)
  # for EMSCRIPTEN GLEW and GLFW are built-in
  include_directories(${PROJECT_SOURCE_DIR}/external/pmp-library/external/glfw/include)
  include_directories(${PROJECT_SOURCE_DIR}/external/pmp-library/external/glew/include)
endif()


# set default compiler flags
if((UNIX OR APPLE) AND CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pedantic -Wall -Wextra")
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" AND NOT EMSCRIPTEN)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pedantic -Wall -Wextra -Wno-unused-parameter")
endif()

if(WIN32)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_USE_MATH_DEFINES -DNOMINMAX -D_CRT_SECURE_NO_WARNINGS")
endif()

if (EMSCRIPTEN)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -s USE_GLFW=3 -s WASM=1 --no-heap-copy -s USE_WEBGL2=1 -s ALLOW_MEMORY_GROWTH=1 -s MINIFY_HTML=0 -s DISABLE_EXCEPTION_CATCHING=0")
    set(CMAKE_EXECUTABLE_SUFFIX ".html")
endif()

# make a release build by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()


if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
endif()


# which directories to process
add_subdirectory(src/)
