
set( gtest_geosx_tests
     testMimeticInnerProducts.cpp )

set( dependencyList mainInterface gtest ${parallelDeps} )

foreach(test ${gtest_geosx_tests})
  get_filename_component( test_name ${test} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )

  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.3.1)
      target_compile_options(${test_name} PRIVATE "-Wno-stringop-overflow")
    endif()
  endif()

  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )
endforeach()
