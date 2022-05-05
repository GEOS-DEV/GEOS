#
# Specify list of tests
#


# the following two lists must have equal length:
# they specify the couples (python script, xml file) corresponding to each test
set( gtest_pygeosx_tests
	testPygeosxAcoustic.py
   )

set( gtest_pygeosx_xmls
        acous_5x5x5.xml
   )

set( dependencyList pygeosx gtest )

if ( GEOSX_BUILD_SHARED_LIBS )
  set (dependencyList ${dependencyList} geosx_core )
else()
  set (dependencyList ${dependencyList} ${geosx_core_libs} )
endif()

if( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

list( LENGTH gtest_pygeosx_tests ntestsp1 )
math( EXPR ntests "${ntestsp1} - 1" )

foreach( itest RANGE ${ntests} )
    list( GET gtest_pygeosx_tests ${itest} pytest )
    list( GET gtest_pygeosx_xmls ${itest} xml )
    get_filename_component( test_script_name ${pytest} NAME_WE )
    get_filename_component( test_xml_name ${xml} NAME_WE )
    string( CONCAT test_name ${test_script_name} "_" ${test_xml_name} )
    message( ${test_name} )
    blt_add_test( NAME ${test_name}
	    COMMAND ${CMAKE_BINARY_DIR}/lib/PYGEOSX/bin/python ${CMAKE_CURRENT_LIST_DIR}/${pytest} -i ${CMAKE_CURRENT_LIST_DIR}/${xml}
                  )
endforeach()
