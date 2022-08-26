file(GLOB SRCS ./*.cpp)
file(GLOB HDRS ./*.h)

add_library(VolumeMesh ${SRCS} ${HDRS})
target_link_libraries(VolumeMesh OpenVolumeMesh pmp_vis)