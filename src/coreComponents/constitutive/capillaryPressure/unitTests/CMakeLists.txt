# Specify list of capillary pressure tests
set( gtest_cappres_tests
     testCapillaryPressure.cpp )

# Base dependency list (mirrors relativePermeability tests plus functions usage)
set( dependencyList gtest blas lapack constitutive mainInterface ${parallelDeps} )

if( ENABLE_CUDA AND ENABLE_CUDA_NVTOOLSEXT )
  list( APPEND dependencyList CUDA::nvToolsExt )
endif()

foreach( test ${gtest_cappres_tests} )
  get_filename_component( test_name ${test} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )
  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )
endforeach()
