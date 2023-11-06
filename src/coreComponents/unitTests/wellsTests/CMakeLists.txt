#
# Specify list of tests
#

set( gtest_geosx_tests
     testReservoirSinglePhaseMSWells.cpp
     testWellEnums.cpp
   )

set( dependencyList ${parallelDeps} gtest )

if ( GEOSX_BUILD_SHARED_LIBS )
  set (dependencyList ${dependencyList} geosx_core ${parallelDeps} )
else()
  set (dependencyList ${dependencyList} ${geosx_core_libs} ${parallelDeps} )
endif()

if( ENABLE_PVTPackage )
    list( APPEND gtest_geosx_tests
          testReservoirCompositionalMultiphaseMSWells.cpp )

    set( dependencyList ${dependencyList} PVTPackage )
endif()


#
# Add gtest C++ based tests
#
foreach(test ${gtest_geosx_tests})
  get_filename_component( test_name ${test} NAME_WE )

  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )

  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )
endforeach()

# For some reason, BLT is not setting CUDA language for these source files
if ( ENABLE_CUDA )
  set_source_files_properties( ${gtest_geosx_tests} PROPERTIES LANGUAGE CUDA )
endif()
