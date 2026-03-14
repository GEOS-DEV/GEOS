# Specify list of tests
set( gtest_geosx_tests
     testThermalCompMultiphaseFlow.cpp
     testThermalSinglePhaseFlow.cpp
     testTransmissibility.cpp
     testImmiscibleMultiphaseFlow.cpp
     testSinglePhaseMFDPolyhedral.cpp
     testSinglePhaseReactiveTransport.cpp )

if( ENABLE_PVTPackage )
    list( APPEND gtest_geosx_tests
          testCompMultiphaseFlow.cpp
          testCompMultiphaseFlowHybrid.cpp
          testReactiveCompositionalMultiphaseOBL.cpp
          testCompositionalMultiPhaseMFDPolyhedral.cpp )
endif()

set( tplDependencyList ${parallelDeps} gtest )

set( dependencyList mainInterface testingUtilities )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

# Copy mesh files to the test output directory once, via a dedicated target,
# to avoid Ninja race conditions when multiple test targets run POST_BUILD
# copy commands for the same destination files in parallel.
set( fluidflow_mesh_files
     polyhedral_voronoi_regular.vtu
     polyhedral_voronoi_lattice.vtu
     polyhedral_voronoi_complex.vtu
     hex_pyr_tet_nested_mixed.vtu )

set( fluidflow_mesh_outputs )
foreach( mesh_file ${fluidflow_mesh_files} )
    set( src "${CMAKE_CURRENT_SOURCE_DIR}/${mesh_file}" )
    set( dst "${TEST_OUTPUT_DIRECTORY}/${mesh_file}" )
    list( APPEND fluidflow_mesh_outputs ${dst} )
    add_custom_command(
        OUTPUT ${dst}
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${src} ${dst}
        DEPENDS ${src}
        COMMENT "Copying ${mesh_file} to test output directory"
    )
endforeach()

add_custom_target( fluidflow_copy_mesh_files ALL
    DEPENDS ${fluidflow_mesh_outputs} )

# Add gtest C++ based tests
foreach(test ${gtest_geosx_tests})
  get_filename_component( test_name ${test} NAME_WE )

  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${decoratedDependencies} ${tplDependencyList} )

  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )

  set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "TEST_BINARY_DIR=$<TARGET_FILE_DIR:${test_name}>")

  # Add TEST_BINARY_DIR as a compile definition
  target_compile_definitions(${test_name} PRIVATE TEST_BINARY_DIR=\"$<TARGET_FILE_DIR:${test_name}>\")

  # Ensure mesh files are copied before any test executable is linked
  add_dependencies( ${test_name} fluidflow_copy_mesh_files )

endforeach()

# testFixedDimHydrostaticEquilibrium: serial (single-rank) variant is always built.

# --- Serial variant (always built) ---
set( test "testFixedDimHydrostaticEquilibrium.cpp" )
get_filename_component( test_name ${test} NAME_WE )

blt_add_executable( NAME ${test_name}
                    SOURCES ${test}
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${decoratedDependencies} ${tplDependencyList} )

geos_add_test( NAME ${test_name}
               COMMAND ${test_name} )

set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "TEST_BINARY_DIR=$<TARGET_FILE_DIR:${test_name}>")
target_compile_definitions(${test_name} PRIVATE TEST_BINARY_DIR=\"$<TARGET_FILE_DIR:${test_name}>\")

# For some reason, BLT is not setting CUDA language for these source files
if ( ENABLE_CUDA )
  set_source_files_properties( ${gtest_geosx_tests} PROPERTIES LANGUAGE CUDA )
endif()
