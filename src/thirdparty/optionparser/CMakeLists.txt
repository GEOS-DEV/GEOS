
message("adding optionparser library")
blt_register_library( 
    NAME optionparser
    INCLUDES optionparser/optionparser.h
    )

set(optionparserheaders optionparser/optionparser.h)
set(optionparserheaderdest ${HEADER_INCLUDES_DIRECTORY}/optionparser)
blt_copy_headers_target(NAME        optionparser
                        HEADERS     ${optionparserheaders}
                        DESTINATION ${optionparserheaderdest})


# The tests throw errors on compilation
# set(optionparserTests
#     tests/testOptionParser.cpp
#     tests/testparse.cpp
#    )

# #
# # Add gtest C++ based tests
# #
# foreach(test ${optionparserTests})
#     get_filename_component( test_name ${test} NAME_WE )
#     blt_add_executable( NAME ${test_name}
#                         SOURCES ${test}
#                         OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
#                         SYSTEM_DEPENDS_ON gtest )
                        
#     blt_add_test( NAME ${test_name} 
#                   COMMAND ${test_name}
#                   )
# endforeach()
