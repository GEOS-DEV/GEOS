# Specify list of tests
set( mesh_tests
     testCellElementRegionSelector.cpp
     testMeshAdjacency.cpp
     testMeshObjectPath.cpp
     testComputationalGeometry.cpp
     testGeometricObjects.cpp
     testGraphColoring.cpp
     testMeshEnums.cpp )
     

set( mesh_tests_mpi
     testGraphColoringMPI.cpp )



set( dependencyList blas lapack gtest mesh ${parallelDeps} )

if( ENABLE_VTK )
  blt_add_executable( NAME testVTKRefinementTemplates
                      SOURCES testVTKRefinementTemplates.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON gtest vtkRefinementTemplates VTK::IOLegacy VTK::IOXML )
  target_compile_definitions( testVTKRefinementTemplates PRIVATE
    VTK_REFINEMENT_FIXTURE_DIR="${PROJECT_SOURCE_DIR}/coreComponents/integrationTests/meshTests" )
  geos_add_test( NAME testVTKRefinementTemplates COMMAND testVTKRefinementTemplates )
  blt_add_executable( NAME testVTKRefinementFields
                      SOURCES testVTKRefinementFields.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON gtest vtkRefinementFields )
  geos_add_test( NAME testVTKRefinementFields COMMAND testVTKRefinementFields )
  blt_add_executable( NAME testVTKRefinementCommunication
                      SOURCES testVTKRefinementCommunication.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON gtest vtkRefinementCommunication vtkRefinementFields ${parallelDeps} )
  target_include_directories( testVTKRefinementCommunication PRIVATE ${PROJECT_SOURCE_DIR}/coreComponents )
  geos_add_test( NAME testVTKRefinementCommunication COMMAND testVTKRefinementCommunication )
  blt_add_executable( NAME testVTKUniformRefinement
                      SOURCES testVTKUniformRefinement.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON gtest vtkUniformRefinement VTK::IOLegacy ${parallelDeps} )
  target_compile_definitions( testVTKUniformRefinement PRIVATE
    VTK_REFINEMENT_FIXTURE_DIR="${PROJECT_SOURCE_DIR}/coreComponents/integrationTests/meshTests" )
  target_include_directories( testVTKUniformRefinement PRIVATE ${PROJECT_SOURCE_DIR}/coreComponents )
  geos_add_test( NAME testVTKUniformRefinement COMMAND testVTKUniformRefinement )
  set_tests_properties( testVTKUniformRefinement PROPERTIES TIMEOUT 60 )
  if( ENABLE_MPI )
    foreach( ranks 2 4 8 )
      geos_add_test( NAME testVTKRefinementCommunication_${ranks}ranks
                     COMMAND testVTKRefinementCommunication NUM_MPI_TASKS ${ranks} )
      set_tests_properties( testVTKRefinementCommunication_${ranks}ranks PROPERTIES TIMEOUT 60 )
      geos_add_test( NAME testVTKUniformRefinement_${ranks}ranks
                     COMMAND testVTKUniformRefinement NUM_MPI_TASKS ${ranks} )
      set_tests_properties( testVTKUniformRefinement_${ranks}ranks PROPERTIES TIMEOUT 60 )
    endforeach()
  endif()
endif()

# Add gtest C++ based tests
foreach(test ${mesh_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} )

    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} )
endforeach()






if( ENABLE_MPI )
set( dependencyList blas lapack gtest mesh ${parallelDeps} )

if (ENABLE_TRILINOS)
   list( APPEND dependencyList trilinos )
endif()

set( nranks 4 )

  foreach( test ${mesh_tests_mpi} )
    get_filename_component( file_we ${test} NAME_WE )
    set( test_name ${file_we}_mpi )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} ${tplDependencyList} )

    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} -x ${nranks}
                   NUM_MPI_TASKS ${nranks} )
    # Older g++ version need a hint to use std::filesystem
    if( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" )
        if( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9 )
            target_link_libraries( ${test_name} PRIVATE "stdc++fs" )
        endif ()
    endif ()
    if ( ENABLE_CUDA )
      target_link_libraries( ${test_name} PRIVATE "stdc++fs" )
    endif ()
  endforeach()

  # VTKMeshScattering test — 4 MPI ranks
  if( ENABLE_VTK )
  blt_add_executable( NAME testVTKMeshScattering_mpi
                      SOURCES testVTKMeshScattering.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} ${tplDependencyList} )

  geos_add_test( NAME testVTKMeshScattering_mpi
                 COMMAND testVTKMeshScattering_mpi
                 NUM_MPI_TASKS 4 )
  set_tests_properties( testVTKMeshScattering_mpi PROPERTIES TIMEOUT 60 )
  endif()
endif()
