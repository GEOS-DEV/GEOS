# Unit tests for elementFormulations (migrated from integrationTests)

set( testSources
     testConformingVirtualElementOrder1.cpp )

set( dependencyList mainInterface gtest ${parallelDeps} )

foreach(test ${testSources})
  get_filename_component( test_name ${test} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )

  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )

  if( ( NOT ${ENABLE_CUDA} ) AND ( CMAKE_CXX_COMPILER_ID STREQUAL "Clang" ) )
    set_target_properties(${test_name} PROPERTIES COMPILE_FLAGS "-Wno-everything")
  endif()
endforeach()
