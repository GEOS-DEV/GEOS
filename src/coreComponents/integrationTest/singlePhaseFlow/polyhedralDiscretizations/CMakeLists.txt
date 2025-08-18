cmake_minimum_required(VERSION 3.10)

# Specify list of tests
set( polyhedral_tests
     testPolyhedralDiscretizations.cpp
   )

set( tplDependencyList gtest )
set( dependencyList mainInterface )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                     DEPENDENCIES ${dependencyList} )

# Add GTest C++ based tests
foreach(test ${polyhedral_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${decoratedDependencies} ${tplDependencyList} )
    geos_add_test( NAME ${test_name} COMMAND ${test_name} )
endforeach()

# For CUDA support
if ( ENABLE_CUDA )
  set_source_files_properties( ${polyhedral_tests} PROPERTIES LANGUAGE CUDA )
endif()
