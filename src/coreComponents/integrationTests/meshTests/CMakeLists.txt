# Specify list of tests
set( gtest_geosx_tests
     testMeshGeneration.cpp
     testNeighborCommunicator.cpp
     testElementRegions.cpp
     testSurfaceGenerator.cpp
     testSplitMesh.cpp )

set( gtest_geosx_mpi_tests
     testNeighborCommunicator.cpp )

if( ENABLE_VTK )
  list( APPEND gtest_geosx_tests
        testVTKImport.cpp )
  list( APPEND gtest_geosx_mpi_tests
        testVTKImport.cpp )
endif()

if( ENABLE_VTK )
  set( TEST_MESH_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/meshDirName.hpp.in ${CMAKE_BINARY_DIR}/include/tests/meshDirName.hpp )
endif()

set( tplDependencyList ${parallelDeps} gtest )

set( dependencyList mainInterface )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

# ---------------------------------------------------------------------------
# Mesh files needed by testSurfaceGenerator
# ---------------------------------------------------------------------------
# DFN-market meshes: one file per fracture combination, covering flat / wavy
# and no-boundary-cut / full-span variants for both hex and tet grids.
set( dfn_market_mesh_files
     # ---- flat, no-boundary-cut, hex ----
     fractured_mesh_hex_DFN_1.vtu
     fractured_mesh_hex_DFN_2.vtu
     fractured_mesh_hex_DFN_3.vtu
     fractured_mesh_hex_DFN_12.vtu
     fractured_mesh_hex_DFN_13.vtu
     fractured_mesh_hex_DFN_23.vtu
     fractured_mesh_hex_DFN_123.vtu
     # ---- flat, no-boundary-cut, tet ----
     fractured_mesh_tet_DFN_1.vtu
     fractured_mesh_tet_DFN_2.vtu
     fractured_mesh_tet_DFN_3.vtu
     fractured_mesh_tet_DFN_12.vtu
     fractured_mesh_tet_DFN_13.vtu
     fractured_mesh_tet_DFN_23.vtu
     fractured_mesh_tet_DFN_123.vtu
     # ---- flat, full-span, hex ----
     fractured_full_span_mesh_hex_DFN_1.vtu
     fractured_full_span_mesh_hex_DFN_2.vtu
     fractured_full_span_mesh_hex_DFN_3.vtu
     fractured_full_span_mesh_hex_DFN_12.vtu
     fractured_full_span_mesh_hex_DFN_13.vtu
     fractured_full_span_mesh_hex_DFN_23.vtu
     fractured_full_span_mesh_hex_DFN_123.vtu
     # ---- flat, full-span, tet ----
     fractured_full_span_mesh_tet_DFN_1.vtu
     fractured_full_span_mesh_tet_DFN_2.vtu
     fractured_full_span_mesh_tet_DFN_3.vtu
     fractured_full_span_mesh_tet_DFN_12.vtu
     fractured_full_span_mesh_tet_DFN_13.vtu
     fractured_full_span_mesh_tet_DFN_23.vtu
     fractured_full_span_mesh_tet_DFN_123.vtu
     # ---- wavy, no-boundary-cut, hex ----
     fractured_wavy_mesh_hex_DFN_1.vtu
     fractured_wavy_mesh_hex_DFN_2.vtu
     fractured_wavy_mesh_hex_DFN_3.vtu
     fractured_wavy_mesh_hex_DFN_12.vtu
     fractured_wavy_mesh_hex_DFN_13.vtu
     fractured_wavy_mesh_hex_DFN_23.vtu
     fractured_wavy_mesh_hex_DFN_123.vtu
     # ---- wavy, no-boundary-cut, tet ----
     fractured_wavy_mesh_tet_DFN_1.vtu
     fractured_wavy_mesh_tet_DFN_2.vtu
     fractured_wavy_mesh_tet_DFN_3.vtu
     fractured_wavy_mesh_tet_DFN_12.vtu
     fractured_wavy_mesh_tet_DFN_13.vtu
     fractured_wavy_mesh_tet_DFN_23.vtu
     fractured_wavy_mesh_tet_DFN_123.vtu
     # ---- wavy, full-span, hex ----
     fractured_wavy_full_span_mesh_hex_DFN_1.vtu
     fractured_wavy_full_span_mesh_hex_DFN_2.vtu
     fractured_wavy_full_span_mesh_hex_DFN_3.vtu
     fractured_wavy_full_span_mesh_hex_DFN_12.vtu
     fractured_wavy_full_span_mesh_hex_DFN_13.vtu
     fractured_wavy_full_span_mesh_hex_DFN_23.vtu
     fractured_wavy_full_span_mesh_hex_DFN_123.vtu
     # ---- wavy, full-span, tet ----
     fractured_wavy_full_span_mesh_tet_DFN_1.vtu
     fractured_wavy_full_span_mesh_tet_DFN_2.vtu
     fractured_wavy_full_span_mesh_tet_DFN_3.vtu
     fractured_wavy_full_span_mesh_tet_DFN_12.vtu
     fractured_wavy_full_span_mesh_tet_DFN_13.vtu
     fractured_wavy_full_span_mesh_tet_DFN_23.vtu
     fractured_wavy_full_span_mesh_tet_DFN_123.vtu
     # ---- T-shaped junction, boundary-cutting, hex & tet ----
     t_shaped_wavy_mesh_hex_DFN_t1t2.vtu
     t_shaped_wavy_mesh_tet_DFN_t1t2.vtu
     # ---- Y-shaped junction, boundary-cutting, hex & tet ----
     y_shaped_wavy_mesh_hex_DFN_y1y2y3.vtu
     y_shaped_wavy_mesh_tet_DFN_y1y2y3.vtu
     # ---- miscellaneous, no-boundary-cut, 5 fractures, hex & tet ----
     DFN_5_fractures_hex_binarized.vtu
     DFN_5_fractures_tet_binarized.vtu
)
set( DFN_MARKET_MESH_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dfn_market )

# ---------------------------------------------------------------------------
# Mesh files needed by testSplitMesh (VTK multiblock .vtm with faceBlocks)
# ---------------------------------------------------------------------------
set( split_mesh_files
     DFN_5_fractures_hex_binarized.vtm
     DFN_5_fractures_hex_binarized_domain.vtu
     DFN_5_fractures_hex_binarized_fault.vtu
)
set( SPLIT_MESH_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dfn_market_split_meshes )


# Add gtest C++ based tests
foreach( test ${gtest_geosx_tests} )
  get_filename_component( test_name ${test} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${decoratedDependencies} ${tplDependencyList} )

  # Add compile definition for TEST_BINARY_DIR.
  # The escaped quotes ensure the path is always a quoted string literal
  # regardless of the CMake generator (Makefiles, Ninja, Xcode, etc.).
  target_compile_definitions( ${test_name} PRIVATE TEST_BINARY_DIR=\"$<TARGET_FILE_DIR:${test_name}>\" )

  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )

  set_tests_properties( ${test_name} PROPERTIES
                        ENVIRONMENT "TEST_BINARY_DIR=$<TARGET_FILE_DIR:${test_name}>"
                        LABELS "${test_name}" )
endforeach()

# ---------------------------------------------------------------------------
# Copy DFN-market mesh files into the test output directory at configure time.
# These are static assets committed to the repository, so a one-time copy is
# sufficient and avoids the race condition that arises with Ninja when
# multiple POST_BUILD commands from different targets write the same files
# to the same directory in parallel.
# ---------------------------------------------------------------------------
file( MAKE_DIRECTORY ${TEST_OUTPUT_DIRECTORY} )
foreach( mesh_file ${dfn_market_mesh_files} )
  file( COPY ${DFN_MARKET_MESH_SOURCE_DIR}/${mesh_file}
        DESTINATION ${TEST_OUTPUT_DIRECTORY} )
endforeach()

# Copy split-mesh (.vtm) files if they exist in the source directory.
foreach( mesh_file ${split_mesh_files} )
  if( EXISTS ${SPLIT_MESH_SOURCE_DIR}/${mesh_file} )
    file( COPY ${SPLIT_MESH_SOURCE_DIR}/${mesh_file}
          DESTINATION ${TEST_OUTPUT_DIRECTORY} )
  endif()
endforeach()

if( ENABLE_MPI AND NOT ENABLE_CUDA AND NOT ENABLE_HIP )

  set( nranks 2 )

  foreach( test ${gtest_geosx_mpi_tests} )
    get_filename_component( file_we ${test} NAME_WE )
    set( test_name ${file_we}_mpi )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${decoratedDependencies} ${tplDependencyList} )

    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} -x ${nranks}
                   NUM_MPI_TASKS ${nranks} )
  endforeach()

  # ---------------------------------------------------------------------------
  # testSurfaceGenerator_mpi — SurfaceGeneratorMpiTest suite (4 MPI ranks)
  #
  # Now built from its own source file (testSurfaceGenerator_mpi.cpp).
  # ---------------------------------------------------------------------------
  set( sg_mpi_nranks 4 )
  set( sg_mpi_test_name testSurfaceGenerator_mpi )

  blt_add_executable( NAME       ${sg_mpi_test_name}
                      SOURCES    testSurfaceGenerator_mpi.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${decoratedDependencies} ${tplDependencyList} )

  target_compile_definitions( ${sg_mpi_test_name}
                              PRIVATE TEST_BINARY_DIR=\"$<TARGET_FILE_DIR:${sg_mpi_test_name}>\" )

  geos_add_test( NAME          ${sg_mpi_test_name}
                 COMMAND       ${sg_mpi_test_name}
                 NUM_MPI_TASKS ${sg_mpi_nranks} )

  set_tests_properties( ${sg_mpi_test_name}
                        PROPERTIES
                        ENVIRONMENT "TEST_BINARY_DIR=$<TARGET_FILE_DIR:${sg_mpi_test_name}>"
                        LABELS "${sg_mpi_test_name}" )

  # DFN-market mesh files are already copied at configure time (see above).

  # ---------------------------------------------------------------------------
  # testSplitMesh_mpi — SplitMeshMpiTest suite (4 MPI ranks)
  # ---------------------------------------------------------------------------
  set( sm_mpi_nranks 4 )
  set( sm_mpi_test_name testSplitMesh_mpi )

  blt_add_executable( NAME       ${sm_mpi_test_name}
                      SOURCES    testSplitMesh_mpi.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${decoratedDependencies} ${tplDependencyList} )

  target_compile_definitions( ${sm_mpi_test_name}
                              PRIVATE TEST_BINARY_DIR=\"$<TARGET_FILE_DIR:${sm_mpi_test_name}>\" )

  geos_add_test( NAME          ${sm_mpi_test_name}
                 COMMAND       ${sm_mpi_test_name}
                 NUM_MPI_TASKS ${sm_mpi_nranks} )

  set_tests_properties( ${sm_mpi_test_name}
                        PROPERTIES
                        ENVIRONMENT "TEST_BINARY_DIR=$<TARGET_FILE_DIR:${sm_mpi_test_name}>"
                        LABELS "${sm_mpi_test_name}" )

endif()
