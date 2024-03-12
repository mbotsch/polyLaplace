cmake_minimum_required(VERSION 3.0.2)
project(PolygonLaplacians)

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

# setup clang format
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# compile PMP library
set(PMP_BUILD_APPS     OFF CACHE BOOL "")
set(PMP_BUILD_EXAMPLES OFF CACHE BOOL "")
set(PMP_BUILD_TESTS    OFF CACHE BOOL "")
set(PMP_BUILD_DOCS     OFF CACHE BOOL "")


# add include directories
include_directories(${PROJECT_SOURCE_DIR}/external/pmp-library/src)
include_directories(${PROJECT_SOURCE_DIR}/external/pmp-library/external/imgui)
include_directories(${PROJECT_SOURCE_DIR}/external/pmp-library/external/eigen)
include_directories(${PROJECT_SOURCE_DIR}/external/spectra/include)
include_directories(${PROJECT_SOURCE_DIR}/external/OpenVolumeMesh/src)
include_directories(${PROJECT_SOURCE_DIR}/external/TinyAD/include)

add_definitions(-DINCLUDE_TEMPLATES)

if (NOT EMSCRIPTEN)
  # for EMSCRIPTEN GLEW and GLFW are built-in
  include_directories(${PROJECT_SOURCE_DIR}/external/pmp-library/external/glfw/include)
  include_directories(${PROJECT_SOURCE_DIR}/external/pmp-library/external/glew/include)
endif()


if(WIN32)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_USE_MATH_DEFINES -DNOMINMAX -D_CRT_SECURE_NO_WARNINGS -fopenmp")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

if (EMSCRIPTEN)
  add_compile_options(--no-heap-copy)
  add_link_options("SHELL:-sWASM=1 -sUSE_WEBGL2=1 -sUSE_GLFW=3 -sALLOW_MEMORY_GROWTH=1 -sALLOW_TABLE_GROWTH=1 -sSTACK_SIZE=5MB -sDISABLE_EXCEPTION_CATCHING=0")
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
add_subdirectory(external/TinyAD)
add_subdirectory(external/pmp-library)
add_subdirectory(external/OpenVolumeMesh)
add_subdirectory(src/)

# add clang-format target if program is found
include(clang-format)