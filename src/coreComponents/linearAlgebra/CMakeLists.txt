message( "Entering /src/coreComponents/linearAlgebra/CMakeLists.txt")
#
# Specify all sources
#
set( linearAlgebra_sources
     interfaces/BlasLapackLA.cpp
     solvers/BiCGSTABsolver.cpp
     solvers/BlockPreconditioner.cpp
     solvers/CGsolver.cpp
     solvers/GMRESsolver.cpp
     solvers/KrylovSolver.cpp
     solvers/SeparateComponentPreconditioner.cpp
     utilities/LAIHelperFunctions.cpp
     DofManager.cpp )

if( BUILD_OBJ_LIBS )
  set( dependencyList common blas lapack RAJA managers)
else()
  set( dependencyList common blas lapack RAJA )
endif()


if( ENABLE_MPI )
    list( APPEND dependencyList mpi )
endif()

list( APPEND linearAlgebra_sources
      interfaces/direct/SuiteSparse.cpp
      interfaces/direct/SuperLU_Dist.cpp )

list( APPEND dependencyList suitesparse superlu_dist )

if( ENABLE_TRILINOS )
    list( APPEND linearAlgebra_sources
          interfaces/trilinos/EpetraMatrix.cpp
          interfaces/trilinos/EpetraVector.cpp
          interfaces/trilinos/TrilinosPreconditioner.cpp
          interfaces/trilinos/TrilinosSolver.cpp
          interfaces/trilinos/TrilinosInterface.cpp
          interfaces/trilinos/EpetraSuiteSparse.cpp
          interfaces/trilinos/EpetraSuperLU_Dist.cpp )

    list( APPEND dependencyList trilinos )

endif()

if( ENABLE_HYPRE )
    list( APPEND linearAlgebra_sources
          interfaces/hypre/HypreMatrix.cpp
          interfaces/hypre/HypreVector.cpp
          interfaces/hypre/HyprePreconditioner.cpp
          interfaces/hypre/HypreSolver.cpp
          interfaces/hypre/HypreInterface.cpp
          interfaces/hypre/HypreSuiteSparse.cpp
          interfaces/hypre/HypreSuperLU_Dist.cpp )

    list( APPEND dependencyList hypre )

endif()

if( ENABLE_PETSC )
    list( APPEND linearAlgebra_sources
          interfaces/petsc/PetscVector.cpp
          interfaces/petsc/PetscMatrix.cpp
          interfaces/petsc/PetscPreconditioner.cpp
          interfaces/petsc/PetscSolver.cpp
          interfaces/petsc/PetscInterface.cpp
          interfaces/petsc/PetscSuiteSparse.cpp
          interfaces/petsc/PetscSuperLU_Dist.cpp )

    list( APPEND dependencyList petsc )

endif()

if( ENABLE_OPENMP )
    set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
    set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  linearAlgebra
                 SOURCES               ${linearAlgebra_sources}
                 DEPENDS_ON            ${dependencyList}
               )

target_include_directories( linearAlgebra PUBLIC ${CMAKE_CURRENT_LIST_DIR} )

geosx_add_code_checks( PREFIX linearAlgebra )

add_subdirectory( unitTests )

message( "Leaving /src/coreComponents/linearAlgebra/CMakeLists.txt")
