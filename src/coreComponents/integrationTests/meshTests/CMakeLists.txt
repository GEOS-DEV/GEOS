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
# Mesh files needed by testSurfaceGenerator (dfn_market — single .vtu per mesh)
# ---------------------------------------------------------------------------
set( DFN_MARKET_MESH_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dfn_market )
file( GLOB dfn_market_mesh_files "${DFN_MARKET_MESH_SOURCE_DIR}/*.vtu" )

# ---------------------------------------------------------------------------
# Mesh files needed by testSplitMesh (dfn_market_split_meshes — .vtm + companion .vtu)
# ---------------------------------------------------------------------------
set( SPLIT_MESH_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dfn_market_split_meshes )
file( GLOB split_mesh_files
      "${SPLIT_MESH_SOURCE_DIR}/*.vtm"
      "${SPLIT_MESH_SOURCE_DIR}/*.vtu" )


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
# Copy all mesh files into the test output directory at configure time.
# file(GLOB ...) returns absolute paths so no source-dir prefix is needed.
# ---------------------------------------------------------------------------
file( MAKE_DIRECTORY ${TEST_OUTPUT_DIRECTORY} )

foreach( mesh_file ${dfn_market_mesh_files} )
  file( COPY ${mesh_file}
        DESTINATION ${TEST_OUTPUT_DIRECTORY} )
endforeach()

foreach( mesh_file ${split_mesh_files} )
  file( COPY ${mesh_file}
        DESTINATION ${TEST_OUTPUT_DIRECTORY} )
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
