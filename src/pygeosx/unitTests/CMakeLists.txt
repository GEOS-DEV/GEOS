#
# Specify list of tests
#


# the following two lists must have equal length:
# they specify the couples (python script, xml file) corresponding to each test
set( pythonFiles
	testPygeosxAcoustic.py
   testWaveEquationAdjoint.py
   )

set( xmlFiles
        acous_5x5x5.xml
        twoLayers2000x2000x2000.xml
   )

set( dependencyList pygeosx )

if ( GEOSX_BUILD_SHARED_LIBS )
  set (dependencyList ${dependencyList} geosx_core )
else()
  set (dependencyList ${dependencyList} ${geosx_core_libs} )
endif()

if( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

# This macro parse a test file and extract test_ function from it to pass it to ctest
macro(parse_test_file filename)
  set(filename ${filename})
  file(READ "${filename}" contents)
  string(REPLACE "\n" ";" rows ${contents})
  foreach(line ${rows})
    string(REGEX MATCHALL "^ *def *test_([A-Za-z_0-9])+ *" found ${line})
    if (found)
      string(REGEX REPLACE "^ *def *test_(([A-Za-z_0-9])+) *" "\\1" name "${found}")
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pytest_template.sh.in
        ${CMAKE_CURRENT_BINARY_DIR}/${name}.sh
        FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
      blt_add_test(NAME ${name} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${name}.sh)
    endif()
  endforeach()
endmacro(parse_test_file)

foreach( xmlFile IN LISTS ${xmlFiles})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${xmlFile} ${CMAKE_CURRENT_BINARY_DIR}/${xmlFile} COPYONLY)
endforeach()

foreach( pythonFile ${pythonFiles})
  parse_test_file(${pythonFile})
endforeach()
