#
# Specify list of tests
#

set( dataRepository_tests
     testDefaultValue.cpp
     testObjectCatalog.cpp
     testRestartBasic.cpp
     testRestartExtended.cpp
     testWrapper.cpp
     testXmlWrapper.cpp
     testBufferOps.cpp
     testPacking.cpp
     testWrapperHelpers.cpp
   )

set( dependencyList gtest conduit )

if ( GEOSX_BUILD_SHARED_LIBS )
  set (dependencyList ${dependencyList} geosx_core)
else()
  set (dependencyList ${dependencyList} ${geosx_core_libs} )
endif()

if ( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

if ( ENABLE_CALIPER )
  set( dependencyList ${dependencyList} caliper adiak )
endif()

#
# Add gtest C++ based tests
#
foreach(test ${dataRepository_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList}
                        )

    blt_add_test( NAME ${test_name}
                  COMMAND ${test_name} )
endforeach()

#   blt_add_executable( NAME testWrapper
#                       SOURCES testWrapper.cpp
#                       OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
#                       DEPENDS_ON dataRepository managers gtest conduit
#                       )

#   blt_add_test( NAME testWrapper
#                 COMMAND testWrapper )