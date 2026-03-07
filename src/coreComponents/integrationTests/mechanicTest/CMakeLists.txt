# Specify list of tests
set( tplDependencyList ${parallelDeps} gtest )

set( dependencyList mainInterface testingUtilities )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

# Mesh files needed by all tests (sourced from the canonical dfn_market directory)
set( DFN_MARKET_MESH_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../meshTests/dfn_market )

set( mesh_files
     fractured_mesh_hex_DFN_1.vtu
     fractured_mesh_hex_DFN_12.vtu
     fractured_mesh_hex_DFN_123.vtu
     fractured_mesh_tet_DFN_1.vtu
     fractured_mesh_tet_DFN_123.vtu
     fractured_wavy_mesh_hex_DFN_1.vtu
     fractured_wavy_mesh_hex_DFN_123.vtu
     fractured_wavy_mesh_tet_DFN_1.vtu
     fractured_wavy_mesh_tet_DFN_123.vtu
     fractured_full_span_mesh_hex_DFN_1.vtu
     fractured_full_span_mesh_hex_DFN_123.vtu
     fractured_full_span_mesh_tet_DFN_1.vtu
     fractured_full_span_mesh_tet_DFN_123.vtu
     t_shaped_wavy_mesh_hex_DFN_t1t2.vtu
     t_shaped_wavy_mesh_tet_DFN_t1t2.vtu
     y_shaped_wavy_mesh_hex_DFN_y1y2y3.vtu
     y_shaped_wavy_mesh_tet_DFN_y1y2y3.vtu
     DFN_5_fractures_hex_binarized.vtu
     DFN_5_fractures_tet_binarized.vtu
)

# testFEMConsistency: multi-rank MPI variant (4 ranks) on CPU-only builds,
# serial (single-rank) variant on GPU or CPU builds.
if( ENABLE_MPI AND NOT ENABLE_CUDA AND NOT ENABLE_HIP )
    set( test "testFEMConsistency_mpi.cpp" )
    # The FractureStressCasesPartitioned suite requires 4 ranks (e.g., 1x1x4 partitioning)
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

    foreach( mesh_file ${mesh_files} )
      add_custom_command(
          TARGET ${test_name} POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy_if_different
                  ${DFN_MARKET_MESH_SOURCE_DIR}/${mesh_file}
                  $<TARGET_FILE_DIR:${test_name}>
      )
    endforeach()
else()
    # Serial (single-rank) variant — used on GPU builds or non-MPI builds
    set( test "testFEMConsistency.cpp" )
    get_filename_component( test_name ${test} NAME_WE )

    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${decoratedDependencies} ${tplDependencyList} )

    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} )

    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "TEST_BINARY_DIR=$<TARGET_FILE_DIR:${test_name}>")
    target_compile_definitions(${test_name} PRIVATE TEST_BINARY_DIR=\"$<TARGET_FILE_DIR:${test_name}>\")

    foreach( mesh_file ${mesh_files} )
      add_custom_command(
          TARGET ${test_name} POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy_if_different
                  ${DFN_MARKET_MESH_SOURCE_DIR}/${mesh_file}
                  $<TARGET_FILE_DIR:${test_name}>
      )
    endforeach()
endif()
