# Specify list of tests
set( gtest_geosx_tests
     testWavePropagation.cpp
     testWavePropagationDG.cpp
     testWavePropagationTaper.cpp
     testWavePropagationQ2.cpp
     testWavePropagationElasticFirstOrder.cpp
     testWavePropagationDAS.cpp
     testWavePropagationElasticVTI.cpp
     testWavePropagationElasticTTI.cpp
     testWavePropagationAttenuationElastic.cpp
     testWavePropagationAttenuationAcoustic.cpp
     testWavePropagationAttenuationAcousticVTI.cpp
     testWavePropagationAcousticFirstOrder.cpp
     testWavePropagationAdjoint1.cpp
   )

set( tplDependencyList ${parallelDeps} gtest )

set( dependencyList mainInterface )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

# Add gtest C++ based tests
if( ENABLE_HIP )
  # The wave solvers currently select no HIP finite-element types
  # (WaveSolverTypeDefSEM/DG.hpp), and their setup consequently reaches
  # explicitly unsupported finite-element dispatch paths. Do not register
  # tests that cannot execute on HIP until those element formulations exist.
  message( STATUS "Skipping wave propagation integration tests for HIP: finite-element HIP support is unavailable." )
else()
foreach(test ${gtest_geosx_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${decoratedDependencies} ${tplDependencyList} )

    # Guard to prevent GCC (version 8) from giving warnings due
    # to some sort of possible conversion from int to long unsigned.
    # See also discussion here: https://github.com/GEOS-DEV/LvArray/pull/250
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.3.1)
        target_compile_options(${test_name} PRIVATE "-Wno-stringop-overflow")
      endif()
    endif()

    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} )

endforeach()
endif()

# For some reason, BLT is not setting CUDA language for these source files
if ( ENABLE_CUDA )
  set_source_files_properties( ${gtest_geosx_tests} PROPERTIES LANGUAGE CUDA )
endif()
