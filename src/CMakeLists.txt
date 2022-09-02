file(GLOB SOURCES ./Surface/*.cpp ./HarmonicBase3D2D/*.cpp ./Volume/*.cpp )
file(GLOB HEADERS ./Surface/*.h ./HarmonicBase3D2D/*.hpp ./Volume/*.h )

add_subdirectory(VolumeMesh)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/VolumeMesh)

SET(OBJS
        /home/abunge/programs/libigl/external/triangle/triangle.o
        )

if (EMSCRIPTEN)

    add_library(poly STATIC ${HEADERS} ${SOURCES})
    target_link_libraries(poly pmp_vis VolumeMesh)

    add_executable(run_volume_tests run_volume_tests.cpp)
    target_link_libraries(run_volume_tests ${OBJS} poly VolumeMesh cholmod)

    add_executable(run_surface_tests run_surface_tests.cpp)
    target_link_libraries(run_surface_tests ${OBJS} poly)

    add_executable(surface_viewer surface_viewer_main.cpp  SurfaceViewer.cpp)
    target_link_libraries(surface_viewer ${OBJS} poly triangle)

    add_executable(volume_viewer volume_viewer_main.cpp VolumeViewer.cpp y)
    target_link_libraries(volume_viewer ${OBJS} poly VolumeMesh cholmod)

    set_target_properties(diamond_laplacian PROPERTIES LINK_FLAGS "--shell-file ${PROJECT_SOURCE_DIR}/external/pmp-library/src/apps/data/shell.html")

else ()

    add_library(poly STATIC ${HEADERS} ${SOURCES})
    target_link_libraries(poly pmp_vis VolumeMesh)

    add_executable(run_volume_tests run_volume_tests.cpp)
    target_link_libraries(run_volume_tests ${OBJS} poly VolumeMesh cholmod)

    add_executable(run_surface_tests run_surface_tests.cpp)
    target_link_libraries(run_surface_tests ${OBJS} pmp poly cholmod)

    add_executable(surface_viewer surface_viewer_main.cpp  SurfaceViewer.cpp SurfaceViewer.h )
    target_link_libraries(surface_viewer ${OBJS} pmp pmp_vis poly cholmod )

    add_executable(volume_viewer volume_viewer_main.cpp VolumeViewer.cpp VolumeViewer.h)
    target_link_libraries(volume_viewer ${OBJS} poly VolumeMesh cholmod)
endif ()

