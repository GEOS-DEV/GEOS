#
# Specify list of tests
#
set( LAI_tests
     BlasLapack
     Matrices
     Vectors
     ExternalSolvers
     KrylovSolvers )

set( nranks 2 )

#
# Add gtest C++ based tests
#
set( dependencyList gtest linearAlgebra )

if ( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

foreach(test_name ${LAI_tests})
  set( exec_name test${test_name} )
  blt_add_executable( NAME ${exec_name}
                      SOURCES ${exec_name}.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )
if ( ENABLE_MPI )
  blt_add_test( NAME testLinearAlgebra_${test_name}
                COMMAND ${exec_name} -x ${nranks}
                NUM_MPI_TASKS ${nranks} )
else()
  blt_add_test( NAME testLinearAlgebra_${test_name}
                COMMAND ${exec_name} )
endif()
endforeach()
