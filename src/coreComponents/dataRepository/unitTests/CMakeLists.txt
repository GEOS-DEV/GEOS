# Specify list of tests
set( dataRepository_tests
     testBufferOps.cpp
     testDefaultValue.cpp
     testDeprecatedGroup.cpp
     testDocumentationGeneration.cpp
     testErrorHandling.cpp
     testObjectCatalog.cpp
     testObsoleteGroup.cpp
     testPacking.cpp
     testWrapperHelpers.cpp
     testWrapper.cpp
     testXmlWrapper.cpp )

set( dependencyList mainInterface gtest ${parallelDeps} )

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
