file(GLOB SOURCES ./Surface/*.cpp  ./Volume/*.cpp ./VolumeMesh/*.cpp ./Harmonic/*.cpp )
file(GLOB HEADERS ./Surface/*.h  ./Volume/*.h ./VolumeMesh/*.h ./Harmonic/*.h)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/VolumeMesh)

add_library(poly STATIC ${HEADERS} ${SOURCES})
target_link_libraries(poly OpenVolumeMesh pmp_vis TinyAD)

if (EMSCRIPTEN)

    add_executable(poly_demo poly_demo.cpp)
    target_link_libraries(poly_demo poly)
    set_target_properties(poly_demo PROPERTIES LINK_FLAGS "--shell-file ${PROJECT_SOURCE_DIR}/external/pmp-library/examples/data/shell.html --preload-file ${PROJECT_SOURCE_DIR}/data/dual_Bunny.off@input.off")

else ()

    add_executable(poly_demo poly_demo.cpp)
    target_link_libraries(poly_demo poly)

    add_executable(run_volume_tests run_volume_tests.cpp)
    target_link_libraries(run_volume_tests poly)

    add_executable(run_surface_tests run_surface_tests.cpp)
    target_link_libraries(run_surface_tests poly)

    add_executable(surface_viewer surface_viewer_main.cpp SurfaceViewer.cpp SurfaceViewer.h)
    target_link_libraries(surface_viewer poly )

    add_executable(volume_viewer volume_viewer_main.cpp VolumeViewer.cpp VolumeViewer.h)
    target_link_libraries(volume_viewer poly)

    add_executable(take_timings run_timings.cpp)
    target_link_libraries(take_timings poly)

endif ()

