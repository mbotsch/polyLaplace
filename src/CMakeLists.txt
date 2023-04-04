file(GLOB SOURCES ./Surface/*.cpp  ./Volume/*.cpp ./HarmonicBase3D2D/*.cpp )
file(GLOB HEADERS ./Surface/*.h  ./Volume/*.h ./HarmonicBase3D2D/*.h)

add_subdirectory(VolumeMesh)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/VolumeMesh)


if (EMSCRIPTEN)

    add_library(poly STATIC ${HEADERS} ${SOURCES})
    target_link_libraries(poly pmp_vis VolumeMesh)

    add_executable(run_volume_tests run_volume_tests.cpp)
    target_link_libraries(run_volume_tests poly VolumeMesh)

    add_executable(run_surface_tests run_surface_tests.cpp)
    target_link_libraries(run_surface_tests poly)

    add_executable(surface_viewer surface_viewer_main.cpp  SurfaceViewer.cpp)
    target_link_libraries(surface_viewer poly triangle)

    add_executable(volume_viewer volume_viewer_main.cpp VolumeViewer.cpp)
    target_link_libraries(volume_viewer poly VolumeMesh)

    add_executable(take_timings run_timings.cpp)
    target_link_libraries(take_timings poly VolumeMesh)

else ()

    add_library(poly STATIC ${HEADERS} ${SOURCES})
    target_link_libraries(poly pmp_vis VolumeMesh)

    add_executable(run_volume_tests run_volume_tests.cpp)
    target_link_libraries(run_volume_tests poly VolumeMesh )

    add_executable(run_surface_tests run_surface_tests.cpp)
    target_link_libraries(run_surface_tests pmp poly )

    add_executable(surface_viewer surface_viewer_main.cpp  SurfaceViewer.cpp SurfaceViewer.h )
    target_link_libraries(surface_viewer pmp pmp_vis poly )

    add_executable(volume_viewer volume_viewer_main.cpp VolumeViewer.cpp VolumeViewer.h)
    target_link_libraries(volume_viewer poly VolumeMesh)

    add_executable(take_timings run_timings.cpp)
    target_link_libraries(take_timings poly VolumeMesh)
endif ()

