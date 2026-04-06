# Unit tests for ImmiscibleMultiphase MFD Kernels

set( serial_tests
     testImmiscibleMultiphaseMFDKernels.cpp )

# Add physicsSolversBase explicitly so RTTI (typeinfo) for PhysicsSolverManager is linked.
# fluidFlowSolvers depends on physicsSolversBase, but BLT's DEPENDS_ON is not always transitive for final link.
set( dependencyList mainInterface common gtest )

foreach( test ${serial_tests} )
  get_filename_component( test_name ${test} NAME_WE )
  blt_add_executable( NAME        ${test_name}
                      SOURCES     ${test}
                      OUTPUT_DIR  ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON  ${dependencyList} )
  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )
endforeach()
