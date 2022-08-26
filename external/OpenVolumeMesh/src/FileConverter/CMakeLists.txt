find_package(Boost 1.43.0)
if(Boost_FOUND)

  add_executable(file_converter ovm_converter.cc)
    target_link_libraries(file_converter OpenVolumeMesh::OpenVolumeMesh Boost::boost)

else()
  message(STATUS "Boost not found. Disabling file_converter.")
endif()
