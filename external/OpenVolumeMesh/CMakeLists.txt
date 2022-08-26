cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
cmake_policy(VERSION 3.12.0)


project (OpenVolumeMesh
    VERSION 2.1.0
    LANGUAGES CXX
    HOMEPAGE_URL https://openvolumemesh.org/
    )

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
    set(OVM_STANDALONE_BUILD TRUE)
else()
    set(OVM_STANDALONE_BUILD FALSE)
endif()

set(OVM_ENABLE_APPLICATIONS ${OVM_STANDALONE_BUILD} CACHE BOOL "Build OpenVolumeMesh applications")
set(OVM_ENABLE_UNITTESTS ${OVM_STANDALONE_BUILD} CACHE BOOL "Build OpenVolumeMesh unit tests")
set(OVM_ENABLE_EXAMPLES ${OVM_STANDALONE_BUILD} CACHE BOOL "Build OpenVolumeMesh examples")
set(OVM_BUILD_DOCUMENTATION ${OVM_STANDALONE_BUILD} CACHE BOOL "Build OpenVolumeMesh documentation")


if (OVM_STANDALONE_BUILD)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/bin")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/lib")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/lib")
    set(OVM_TARGET_PREFIX "")
    enable_testing()
    option(STL_DEBUG "Enable STL debug checks. Warning: This must be consistent across the whole build! Not compatible with the libc++ shipped with Apple XCode as of 2019." OFF)
    if (STL_DEBUG)
        # We do not know which STL will be used, defining both is the easist way.

        # https://libcxx.llvm.org/docs/UsingLibcxx.html#libc-configuration-macros
        # https://libcxx.llvm.org/docs/DesignDocs/DebugMode.html#using-debug-mode
        add_definitions(-D_LIBCPP_DEBUG=1)
        # https://gcc.gnu.org/onlinedocs/libstdc%2B%2B/manual/debug_mode_using.html
        add_definitions(-D_GLIBCXX_DEBUG=1)
        add_definitions(-D_GLIBCXX_DEBUG_PEDANTIC=1)
    endif()
    set(OVM_CXX_STANDARD 11 CACHE STRING "C++ standard version to use")
    set_property(CACHE OVM_CXX_STANDARD PROPERTY STRINGS 11 14 17 20)
    set(CMAKE_CXX_STANDARD ${OVM_CXX_STANDARD})
    set(CMAKE_CXX_STANDARD_REQUIRED YES)
    set(CMAKE_CXX_EXTENSIONS OFF)
else()
    # OVM_TARGET_PREFIX is prepended to targets with generic names,
    # i.e., "doc", "unittests"
    set(OVM_TARGET_PREFIX "OVM_")
endif()

set (CMAKE_DEBUG_POSTFIX "d")

add_subdirectory (src)

if (OVM_ENABLE_EXAMPLES)
    add_subdirectory(examples)
endif()

if (OVM_BUILD_DOCUMENTATION)
   add_subdirectory(documentation)
endif()

#export(PACKAGE OpenVolumeMesh)

#list(REMOVE_AT CMAKE_MODULE_PATH -1)
