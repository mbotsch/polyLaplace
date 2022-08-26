SET(SOURCE_FILES
    unittests_smartptr.cc
    unittests_basics.cc
    unittests_files.cc
    unittests_common.cc
    unittests.cc
    unittests_iterators.cc
    unittests_properties.cc)

if (NOT TARGET OpenVolumeMesh::OpenVolumeMesh)
    find_package(OpenVolumeMesh REQUIRED)
endif()

set(TARGET_NAME "${OVM_TARGET_PREFIX}unittests")

enable_testing()
find_package(GTest)

if(GTEST_FOUND)
    include(GoogleTest)

    add_executable(${TARGET_NAME} ${SOURCE_FILES})
    target_link_libraries(${TARGET_NAME}
        OpenVolumeMesh::OpenVolumeMesh
        GTest::GTest
        GTest::Main
        )

    gtest_add_tests(TARGET ${TARGET_NAME}
        EXTRA_ARGS "--gtest_color=yes;--gtest_output=xml"
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/Unittests"
        )

    file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/TestFiles/"
        DESTINATION "${CMAKE_BINARY_DIR}/Unittests"
        )

    set_target_properties(${TARGET_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Unittests
        )

    #[[=
    if ( WIN32 )
        add_custom_command (TARGET ${TARGET_NAME} POST_BUILD
                                COMMAND ${CMAKE_COMMAND} -E
                                    copy_if_different
                                    ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/unittests.exe
                                    ${CMAKE_BINARY_DIR}/Unittests/unittests.exe)
    endif()
    =]]
endif()
