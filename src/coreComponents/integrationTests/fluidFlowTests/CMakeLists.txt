# Specify list of tests
set( gtest_geosx_tests
     testThermalCompMultiphaseFlow.cpp
     testThermalSinglePhaseFlow.cpp
     testTransmissibility.cpp
     testImmiscibleMultiphaseFlow.cpp
     testSinglePhaseMFDPolyhedral.cpp
     testSinglePhaseReactiveTransport.cpp
     testMixedDimSinglePhaseFlow.cpp )

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
    add_custom_command(
        TARGET ${test_name} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_CURRENT_SOURCE_DIR}/../mechanicTest/fractured_mesh_tet_DFN_12.vtu
                $<TARGET_FILE_DIR:${test_name}>
    )

endforeach()

# For some reason, BLT is not setting CUDA language for these source files
if ( ENABLE_CUDA )
  set_source_files_properties( ${gtest_geosx_tests} PROPERTIES LANGUAGE CUDA )
endif()
