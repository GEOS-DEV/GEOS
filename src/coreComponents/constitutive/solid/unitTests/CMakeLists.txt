# Unit tests for constitutive solid models (DamageSpectral etc.)

set( gtest_geosx_tests
     testDamage.cpp 
     testFiniteDifferenceElasticFiniteStrain.cpp )

set( dependencyList gtest mainInterface blas lapack constitutive ${parallelDeps} )

if( ENABLE_CUDA AND ENABLE_CUDA_NVTOOLSEXT )
  list( APPEND dependencyList CUDA::nvToolsExt )
endif()

foreach(test ${gtest_geosx_tests})
  get_filename_component( test_name ${test} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )

  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )
endforeach()
