cmake_minimum_required(VERSION 3.10)

# Specify list of tests
set( polyhedral_tests
     testPolyhedralDiscretizations.cpp
   )

set( tplDependencyList gtest )
set( dependencyList mainInterface )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                     DEPENDENCIES ${dependencyList} )

# Add GTest C++ based tests
foreach(test ${polyhedral_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${decoratedDependencies} ${tplDependencyList} )
    geos_add_test( NAME ${test_name} COMMAND ${test_name} )

    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "TEST_BINARY_DIR=$<TARGET_FILE_DIR:${test_name}>")

    # Add TEST_BINARY_DIR as a compile definition
    target_compile_definitions(${test_name} PRIVATE TEST_BINARY_DIR=\"$<TARGET_FILE_DIR:${test_name}>\")

    # --- Copy mesh file to test binary dir ---
    add_custom_command(
        TARGET ${test_name} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_CURRENT_SOURCE_DIR}/polyhedral_voronoi_regular.vtk
                $<TARGET_FILE_DIR:${test_name}>
    )
    add_custom_command(
        TARGET ${test_name} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_CURRENT_SOURCE_DIR}/polyhedral_voronoi_lattice.vtk
                $<TARGET_FILE_DIR:${test_name}>
    )
    add_custom_command(
        TARGET ${test_name} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_CURRENT_SOURCE_DIR}/polyhedral_voronoi_complex.vtk
                $<TARGET_FILE_DIR:${test_name}>
    )
endforeach()

# For CUDA support
if ( ENABLE_CUDA )
  set_source_files_properties( ${polyhedral_tests} PROPERTIES LANGUAGE CUDA )
endif()
