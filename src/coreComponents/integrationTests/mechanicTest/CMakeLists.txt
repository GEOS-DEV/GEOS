# Specify list of tests
set( tplDependencyList ${parallelDeps} gtest )

set( dependencyList mainInterface testingUtilities )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

# Mesh files needed by all tests (sourced from the canonical dfn_market directory)
set( DFN_MARKET_MESH_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../meshTests/dfn_market )

file( GLOB mesh_files "${DFN_MARKET_MESH_SOURCE_DIR}/*.vtu" )

# ---------------------------------------------------------------------------
# Copy DFN-market mesh files into the test output directory at configure time.
# ---------------------------------------------------------------------------
file( MAKE_DIRECTORY ${TEST_OUTPUT_DIRECTORY} )
foreach( mesh_file ${mesh_files} )
  file( COPY ${mesh_file}
        DESTINATION ${TEST_OUTPUT_DIRECTORY} )
endforeach()

# testFEMConsistency: serial (single-rank) variant is always built.
# Additionally, the multi-rank MPI variant (4 ranks) is built on CPU-only MPI builds.

# --- Serial variant (always built) ---
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

# --- MPI variant (CPU-only MPI builds) ---
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
endif()

# ---------------------------------------------------------------------------
# testFEMConsistencySplitMesh: tests using pre-split (VTM) meshes with
# fractures defined via faceBlocks instead of SurfaceGenerator.
# ---------------------------------------------------------------------------
set( DFN_MARKET_SPLIT_MESH_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../meshTests/dfn_market_split_meshes )

# Glob all VTM and companion VTU files from the split-mesh directory so that
# every mesh available is automatically copied to the test output directory.
file( GLOB split_mesh_files
      "${DFN_MARKET_SPLIT_MESH_SOURCE_DIR}/*.vtm"
      "${DFN_MARKET_SPLIT_MESH_SOURCE_DIR}/*.vtu" )

foreach( mesh_file ${split_mesh_files} )
  file( COPY ${mesh_file}
        DESTINATION ${TEST_OUTPUT_DIRECTORY} )
endforeach()

# --- Serial variant (always built) ---
set( test "testFEMConsistencySplitMesh.cpp" )
get_filename_component( test_name ${test} NAME_WE )

blt_add_executable( NAME ${test_name}
                    SOURCES ${test}
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${decoratedDependencies} ${tplDependencyList} )

geos_add_test( NAME ${test_name}
               COMMAND ${test_name} )

set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "TEST_BINARY_DIR=$<TARGET_FILE_DIR:${test_name}>")
target_compile_definitions(${test_name} PRIVATE TEST_BINARY_DIR=\"$<TARGET_FILE_DIR:${test_name}>\")

# --- MPI variant (CPU-only MPI builds) ---
if( ENABLE_MPI AND NOT ENABLE_CUDA AND NOT ENABLE_HIP )
    set( test "testFEMConsistencySplitMesh_mpi.cpp" )
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
endif()
