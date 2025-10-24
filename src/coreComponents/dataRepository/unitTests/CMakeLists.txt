# Specify list of tests
set( dataRepository_tests
     testDefaultValue.cpp
     testDocumentationGeneration.cpp
     testPacking.cpp
     testWrapper.cpp
     testXmlWrapper.cpp
     testBufferOps.cpp )

set( dependencyList ${parallelDeps} gtest dataRepository )

if( ENABLE_CUDA AND ENABLE_CUDA_NVTOOLSEXT )
  list( APPEND dependencyList CUDA::nvToolsExt )
endif()

# Add gtest C++ based tests
foreach(test ${dataRepository_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} )

    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} )
endforeach()

# Python-based tests
if( ENABLE_PYGEOSX )
    set( pythonTests
         testPyWrapper.cpp )

    foreach(test ${pythonTests})
        get_filename_component( test_name ${test} NAME_WE )

        blt_add_library( NAME ${test_name}
                         SOURCES ${test}
                         DEPENDS_ON pygeosx ${dependencyList}
                         SHARED TRUE
                         CLEAR_PREFIX TRUE )

        target_include_directories( ${test_name} PUBLIC ${CMAKE_CURRENT_LIST_DIR}/../../src )

        configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}Driver.py 
                        ${CMAKE_CURRENT_BINARY_DIR}/${test_name}Driver.py 
                        COPYONLY )

        add_test( NAME ${test_name}
                  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/${test_name}Driver.py )

        set_tests_properties( ${test_name}
                              PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/lib )

    endforeach()
endif()
