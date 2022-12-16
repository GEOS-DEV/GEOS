#add_subdirectory(sphinx)

if(DOXYGEN_FOUND)
    add_subdirectory(doxygen)
endif()


# Python formatting
if ( ENABLE_YAPF )
     set( python_example_sources )
     file( GLOB_RECURSE python_example_sources "*py" )

     # Note: blt throws an error if sources doesn't include a c-file, so include dummy.cpp
     blt_add_code_checks( PREFIX  python_documentation_examples_yapf_style
                          SOURCES ${python_example_sources} ${CMAKE_SOURCE_DIR}/coreComponents/dummy.cpp
                          YAPF_CFG_FILE ${PROJECT_SOURCE_DIR}/yapf.cfg )
endif()

