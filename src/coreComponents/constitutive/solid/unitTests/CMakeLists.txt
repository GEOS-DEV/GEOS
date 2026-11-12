# Unit tests for constitutive solid models (DamageSpectral etc.)

set( gtest_geosx_tests
     testDamage.cpp
     testSolidThermalExpansion.cpp )

set( dependencyList gtest blas lapack constitutive ${parallelDeps} )

if( ENABLE_CUDA AND ENABLE_CUDA_NVTOOLSEXT )
  list( APPEND dependencyList CUDA::nvToolsExt )
endif()

# testSolidThermalExpansion builds coupled solids from XML through the constitutive catalog: link the whole
# constitutive library so that every catalog entry is registered, not only the ones the test references.
geos_decorate_link_dependencies( LIST wholeConstitutiveDependencies
                                 DEPENDENCIES constitutive )

foreach(test ${gtest_geosx_tests})
  get_filename_component( test_name ${test} NAME_WE )
  set( testDependencies ${dependencyList} )
  if( test_name STREQUAL "testSolidThermalExpansion" )
    set( testDependencies gtest blas lapack ${wholeConstitutiveDependencies} ${parallelDeps} )
  endif()
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${testDependencies} )

  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )
endforeach()
