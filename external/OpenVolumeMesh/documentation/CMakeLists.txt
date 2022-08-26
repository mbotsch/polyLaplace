find_package(Doxygen)

if (DOXYGEN_FOUND)
    set(DOXYGEN_PROJECT_LOGO "${CMAKE_CURRENT_SOURCE_DIR}/html/logo_vci.png")
    set(DOXYGEN_IMAGE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/images")
    set(DOXYGEN_WARN_IF_UNDOCUMENTED "FALSE")
    set(DOXYGEN_CITE_BIB_FILES "${CMAKE_CURRENT_SOURCE_DIR}/references.bib")
    set(DOXYGEN_STRIP_FROM_PATH "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/../src")
    set(DOXYGEN_STRIP_FROM_INC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../src")
    set(DOXYGEN_FILE_PATTERNS "*.cc" "*.hh" "*.docu")
    set(DOXYGEN_EXAMPLE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/examples")
    set(DOXYGEN_EXAMPLE_PATTERNS "*.cc" "*.hh" "*.txt" ".ovm")
    set(DOXYGEN_ENABLED_SECTIONS "OPENVOLUMEMESH_INTERNAL_DOC")
    set(DOXYGEN_GENERATE_TREEVIEW "YES")
    set(DOXYGEN_USE_MATHJAX "YES")
    #set(DOXYGEN_EXTRA_PACKAGS "amsmath" "amsfonts" "amssymb" "biblatex")

    set(DOXYGEN_EXCLUDE # "${CMAKE_BINARY_DIR}"
        # "${CMAKE_CURRENT_SOURCE_DIR}/../examples"
                        "${CMAKE_CURRENT_SOURCE_DIR}/../src/Unittests"
                        )

    if (${VCI_PROJECT_DATADIR})
        set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/Build/${VCI_PROJECT_DATADIR}/Doc")
    else()
        set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/Build/share/doc/OpenVolumeMesh")
    endif()

    doxygen_add_docs("${OVM_TARGET_PREFIX}doc"
        "${CMAKE_CURRENT_SOURCE_DIR}/"
        "${CMAKE_CURRENT_SOURCE_DIR}/../src"
        )

endif()
