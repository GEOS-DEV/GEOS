#add_subdirectory(sphinx)

if(DOXYGEN_FOUND)
    add_subdirectory(doxygen)
endif()


# Python formatting
if ( ENABLE_YAPF )
     set( python_example_sources )
     file( GLOB_RECURSE python_example_sources "*py" )

     blt_add_code_checks( PREFIX  python_examples_${arg_PREFIX}
                          SOURCES ${python_example_sources}
                          YAPF_CFG_FILE ${PROJECT_SOURCE_DIR}/yapf.cfg )
endif()

