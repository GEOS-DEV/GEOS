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

    # --- Copy mesh file to test binary dir ---
    add_custom_command(
        TARGET ${test_name} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_CURRENT_SOURCE_DIR}/polyhedral_voronoi_regular.vtu
                $<TARGET_FILE_DIR:${test_name}>
    )
    add_custom_command(
        TARGET ${test_name} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_CURRENT_SOURCE_DIR}/polyhedral_voronoi_lattice.vtu
                $<TARGET_FILE_DIR:${test_name}>
    )
    add_custom_command(
        TARGET ${test_name} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_CURRENT_SOURCE_DIR}/polyhedral_voronoi_complex.vtu
                $<TARGET_FILE_DIR:${test_name}>
    )
    add_custom_command(
        TARGET ${test_name} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_CURRENT_SOURCE_DIR}/hex_pyr_tet_nested_mixed.vtu
                $<TARGET_FILE_DIR:${test_name}>
    )

endforeach()

# Special handling for testMixedDimSinglePhaseFlow which requires 4 MPI tasks
if( ENABLE_MPI )
    set( test "testMixedDimSinglePhaseFlow.cpp" )
    # The MixedDimPartitionedFlowCases suite requires 4 ranks (e.g., 1x1x4 partitioning)
    set( nranks 4 )

    get_filename_component( test_name ${test} NAME_WE )

    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${decoratedDependencies} ${tplDependencyList} )

    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} -x ${nranks}
                   NUM_MPI_TASKS ${nranks} )

    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "TEST_BINARY_DIR=$<TARGET_FILE_DIR:${test_name}>")
    target_compile_definitions(${test_name} PRIVATE TEST_BINARY_DIR=\"$<TARGET_FILE_DIR:${test_name}>\")

    # Copy necessary mesh files for testMixedDimSinglePhaseFlow
    foreach( mesh IN ITEMS "fractured_mesh_hex_DFN_1.vtu"
                           "fractured_mesh_hex_DFN_12.vtu"
                           "fractured_mesh_hex_DFN_123.vtu" )
      add_custom_command(
          TARGET ${test_name} POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy_if_different
                  ${CMAKE_CURRENT_SOURCE_DIR}/../mechanicTest/${mesh}
                  $<TARGET_FILE_DIR:${test_name}>
      )
    endforeach()
endif()

# For some reason, BLT is not setting CUDA language for these source files
if ( ENABLE_CUDA )
  set_source_files_properties( ${gtest_geosx_tests} PROPERTIES LANGUAGE CUDA )
endif()
