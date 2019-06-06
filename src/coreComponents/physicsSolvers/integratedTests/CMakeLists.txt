#
# Specify list of tests
#

set(testPhysicsSolverPackage1
    testSolidMechanicsLagrangianFEM.cpp
    testLaplaceFEM.cpp
   )

#
# Add gtest C++ based tests
#
foreach(test ${testPhysicsSolverPackage1})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON geosx_core  gtest )
                        
    blt_add_test( NAME ${test_name} 
                  COMMAND ${test_name}
                  )
endforeach()
