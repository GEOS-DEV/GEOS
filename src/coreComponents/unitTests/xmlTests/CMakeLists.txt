#
# Specify list of tests
#

set( gtest_geosx_tests
     testXML.cpp
   )

set( gtest_tests_with_input
     testXMLFile.cpp
   )

set( input_files
     basic_input.xml
     include_input.xml
   )

set( multi_files
     multifile_input/fieldspec.xml
     multifile_input/mesh.xml
     multifile_input/outputs.xml
     multifile_input/solver.xml
   )

set( input_data
     multifile_input/boundary_conditions/constraints.xml
     multifile_input/boundary_conditions/velocity.xml
     multifile_input/boundary_conditions/tables/time.txt
     multifile_input/boundary_conditions/tables/velocity.txt
   )
#
# Add gtest C++ based tests
#
set( dependencyList gtest optionparser )

if ( GEOSX_BUILD_SHARED_LIBS )
  set (dependencyList ${dependencyList} geosx_core)
else()
  set (dependencyList ${dependencyList} ${geosx_core_libs} )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

# create test executables
foreach(test ${gtest_geosx_tests} ${gtest_tests_with_input})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} "-lz"
                      )
endforeach()

# add tests that don't require input files
foreach(test ${gtest_geosx_tests})
  get_filename_component( test_name ${test} NAME_WE )
  blt_add_test( NAME ${test_name}
                COMMAND ${test_name} )
endforeach()

foreach( file ${input_files} ${multi_files} ${input_data} )
  configure_file( ${file} ${CMAKE_CURRENT_BINARY_DIR}/${file} COPYONLY )
endforeach()

# add tests that do require input files
foreach( test ${gtest_tests_with_input} )
  get_filename_component( test_name ${test} NAME_WE )

  # test a single input file
  foreach( input ${input_files} )
    get_filename_component( input_name ${input} NAME_WE )
    blt_add_test( NAME ${test_name}_${input_name}
                  COMMAND ${test_name} -i ${CMAKE_CURRENT_BINARY_DIR}/${input} )
  endforeach()

  # test multiple input files
  foreach( file ${multi_files} )
    set( args ${args} -i ${CMAKE_CURRENT_BINARY_DIR}/${file} )
  endforeach()
  blt_add_test( NAME ${test_name}_multiple_input
                COMMAND ${test_name} ${args} )
endforeach()
