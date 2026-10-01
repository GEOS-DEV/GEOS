#[[
Package: fluidFlowSolvers unitTests
Builds unit test for fluid flow solver statistics
#]]

set( fluidFlow_gtest_tests
     testFlowStatistics.cpp
     testComponentDispatch.cpp )

set( dependencyList mainInterface testingUtilities gtest ${parallelDeps} )

foreach( test_source ${fluidFlow_gtest_tests} )
  get_filename_component( test_name ${test_source} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test_source}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList}
                      FOLDER UnitTests )
  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )
  
  if( ENABLE_MPI AND NOT ENABLE_CUDA AND NOT ENABLE_HIP )
    set( nranks 2 )
    geos_add_test( NAME ${test_name}_mpi
                   COMMAND ${test_name} -x ${nranks}
                   NUM_MPI_TASKS ${nranks} )
  endif()
endforeach()

add_test( NAME testComponentOptions
          COMMAND ${CMAKE_COMMAND} -P ${PROJECT_SOURCE_DIR}/cmake/tests/ComponentOptionsTest.cmake )

if( GEOS_ENABLE_MULTIPHYSICS )
  set( componentSmokeInput "${PROJECT_SOURCE_DIR}/../inputFiles/compositionalMultiphaseWell/compositional_multiphase_wells_9comp_smoke.xml" )
  set( componentSmokeOutput "${CMAKE_CURRENT_BINARY_DIR}/componentLimitSmoke" )
  file( MAKE_DIRECTORY "${componentSmokeOutput}" )
  if( GEOS_MAX_COMPONENTS EQUAL 9 )
    add_test( NAME testNineComponentWells
              COMMAND $<TARGET_FILE:geosx> -i ${componentSmokeInput} -o ${componentSmokeOutput} -n nineComponentWells )
    set_tests_properties( testNineComponentWells PROPERTIES TIMEOUT 120 )
    set( componentDirichletInput "${PROJECT_SOURCE_DIR}/../inputFiles/compositionalMultiphaseWell/compositional_multiphase_wells_9comp_dirichlet_smoke.xml" )
    set( componentDirichletOutput "${CMAKE_CURRENT_BINARY_DIR}/componentDirichletSmoke" )
    file( MAKE_DIRECTORY "${componentDirichletOutput}" )
    add_test( NAME testNineComponentDirichlet
              COMMAND $<TARGET_FILE:geosx> -i ${componentDirichletInput} -o ${componentDirichletOutput} -n nineComponentDirichlet )
    set_tests_properties( testNineComponentDirichlet PROPERTIES TIMEOUT 120 )
  else()
    add_test( NAME testComponentLimitRejection
              COMMAND ${CMAKE_COMMAND} -DGEOS_EXECUTABLE=$<TARGET_FILE:geosx>
                      -DINPUT_FILE=${componentSmokeInput} -DOUTPUT_DIRECTORY=${componentSmokeOutput}
                      -P ${PROJECT_SOURCE_DIR}/cmake/tests/ComponentLimitRejectionTest.cmake )
    set_tests_properties( testComponentLimitRejection PROPERTIES TIMEOUT 120 )
  endif()
endif()
