#
# Specify list of tests
#

set(gtest_geosx_test testSinglePhaseFlow.cpp)

get_filename_component(exec_name ${gtest_geosx_test} NAME_WE)
message("processing test_name ${exec_name}")
blt_add_executable(NAME ${exec_name}
                   SOURCES ${gtest_geosx_test}
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${geosx_core_libs} ${thirdPartyLibs} gtest
                   )

set(gtest_geosx_inputs
    2D_10x10_incompr_cos.xml
    2D_100x100_incompr_cos.xml
    )

foreach(input ${gtest_geosx_inputs})
    get_filename_component(test_name ${input} NAME_WE)
    blt_add_test(NAME ${exec_name}/${test_name}
                 COMMAND ${exec_name} -i ${CMAKE_CURRENT_LIST_DIR}/${input}
            )
endforeach()