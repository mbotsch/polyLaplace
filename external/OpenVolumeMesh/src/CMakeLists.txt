

SET(SOURCE_FILES
    OpenVolumeMesh/Attribs/InterfaceAttrib.cc
    OpenVolumeMesh/Attribs/OpenVolumeMeshStatus.cc
    OpenVolumeMesh/Attribs/StatusAttrib.cc
    OpenVolumeMesh/Geometry/Vector.cc
    OpenVolumeMesh/FileManager/FileManager.cc
    OpenVolumeMesh/Core/PropertyDefines.cc
    OpenVolumeMesh/Core/OpenVolumeMeshHandle.cc
    OpenVolumeMesh/Core/ResourceManager.cc
    OpenVolumeMesh/Core/BaseEntities.cc
    OpenVolumeMesh/Core/BaseProperty.cc
    OpenVolumeMesh/Core/TopologyKernel.cc
    OpenVolumeMesh/Core/Serializers.cc
    OpenVolumeMesh/Core/Iterators.cc
    OpenVolumeMesh/Core/OpenVolumeMeshBaseProperty.cc
    OpenVolumeMesh/Mesh/TetrahedralMeshIterators.cc
    OpenVolumeMesh/Mesh/HexahedralMeshIterators.cc
    OpenVolumeMesh/Mesh/TetrahedralMeshTopologyKernel.cc
    OpenVolumeMesh/Mesh/HexahedralMeshTopologyKernel.cc
)

set(OVM_ENABLE_DEPRECATED_APIS OFF CACHE BOOL "Support deprecated OpenVolumeMesh APIs")

configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/OpenVolumeMesh/Config/DeprecationConfig.hh.in"
    "${CMAKE_CURRENT_BINARY_DIR}/OpenVolumeMesh/Config/DeprecationConfig.hh"
)

configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/OpenVolumeMesh/Config/Version.hh.in"
    "${CMAKE_CURRENT_BINARY_DIR}/OpenVolumeMesh/Config/Version.hh"
)


# If we are in a subdirectory of another project that has a VCICommon module
# (e.g. OpenFlipper), we will use vci_add_library; otherwise behave "normally".
include(VCICommon OPTIONAL RESULT_VARIABLE IN_VCI_BUILD)
if(IN_VCI_BUILD)
    message("OVM is being built as part of an VCI project, using vci_add_library!")
    if (WIN32 OR EMSCRIPTEN)
        vci_add_library(OpenVolumeMesh STATIC ${SOURCE_FILES})
    else()
        vci_add_library(OpenVolumeMesh SHARED ${SOURCE_FILES})
    endif()
else(IN_VCI_BUILD)
    if (WIN32)
        # Windows .dll build is not supported yet.
        # TODO: turn this on once we have OVM_EXPORT macros in all required places,
        # cf CXX_VISIBILITY_PRESET:
        add_library(OpenVolumeMesh STATIC ${SOURCE_FILES})
    else()
        option(BUILD_SHARED_LIBS "Build shared library (*.dll, *.so, *.dylib) instead of static library (*.a, *.lib)" ON)
        add_library(OpenVolumeMesh ${SOURCE_FILES})
    endif()
endif(IN_VCI_BUILD)

add_library(OpenVolumeMesh::OpenVolumeMesh ALIAS OpenVolumeMesh)

include(GenerateExportHeader)
generate_export_header(OpenVolumeMesh
    BASE_NAME OVM
    EXPORT_FILE_NAME "OpenVolumeMesh/Config/Export.hh"
    DEPRECATED_MACRO_NAME "CMAKE_OVM_DEPRECATED"
    NO_DEPRECATED_MACRO_NAME "CMAKE_OVM_NO_DEPRECATED"
)

set_target_properties (OpenVolumeMesh PROPERTIES
    VERSION ${OpenVolumeMesh_VERSION_MAJOR}.${OpenVolumeMesh_VERSION_MINOR}
    SOVERSION ${OpenVolumeMesh_VERSION_MAJOR}.${OpenVolumeMesh_VERSION_MINOR}
    POSITION_INDEPENDENT_CODE ON
    # TODO: turn this on once we have OVM_EXPORT macros in all required places:
     #CXX_VISIBILITY_PRESET hidden
    )


target_compile_options(OpenVolumeMesh PUBLIC
    "$<$<CXX_COMPILER_ID:MSVC>:/bigobj>"
    "$<$<CXX_COMPILER_ID:MSVC>:/Zc:__cplusplus>" # Enable updated __cplusplus macro
    PRIVATE
    "$<$<CXX_COMPILER_ID:GNU>:-Wall>"
    "$<$<CXX_COMPILER_ID:Clang>:-Wall;-Wno-c++98-compat;-Wno-c++98-compat-pedantic;-Wno-weak-vtables;-Wno-padded;-Wno-sign-conversion;-Wno-missing-prototypes>"
    "$<$<CXX_COMPILER_ID:MSVC>:/W1>"
    # disable MSVC inline warnings that are enabled with /Ob2, which is in the default cmake release cxxflags:
    "$<$<CXX_COMPILER_ID:MSVC>:/wd4710>"
    "$<$<CXX_COMPILER_ID:MSVC>:/wd4711>"
    )

target_compile_definitions(OpenVolumeMesh PUBLIC "INCLUDE_TEMPLATES")

target_include_directories(OpenVolumeMesh
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)


set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/OpenVolumeMesh)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/OpenVolumeMeshConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/OpenVolumeMeshConfig.cmake"
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
    )

write_basic_package_version_file(
    OpenVolumeMeshConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMinorVersion
    )

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/OpenVolumeMeshConfigVersion.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/OpenVolumeMeshConfig.cmake"
    DESTINATION ${INSTALL_CONFIGDIR})


install(DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/OpenVolumeMesh/Config"
    DESTINATION include/OpenVolumeMesh)

# Install Header Files
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/OpenVolumeMesh
    DESTINATION include
    FILES_MATCHING PATTERN "*.hh"
    )

install(TARGETS OpenVolumeMesh
    EXPORT OpenVolumeMeshTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(EXPORT OpenVolumeMeshTargets
    FILE OpenVolumeMeshTargets.cmake
    NAMESPACE OpenVolumeMesh::
    DESTINATION ${INSTALL_CONFIGDIR}
)

export(EXPORT OpenVolumeMeshTargets
    NAMESPACE OpenVolumeMesh::)
#export(PACKAGE OpenVolumeMesh)

if (OVM_ENABLE_APPLICATIONS)
    add_subdirectory(FileConverter)
endif()

if (OVM_ENABLE_UNITTESTS)
    add_subdirectory(Unittests)
endif()
