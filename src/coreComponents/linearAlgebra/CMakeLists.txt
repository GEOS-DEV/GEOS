message( "Entering /src/coreComponents/linearAlgebra/CMakeLists.txt")

#
# Specify all headers
#
set( linearAlgebra_headers
     common.hpp
     interfaces/BlasLapackFunctions.h
     interfaces/BlasLapackLA.hpp
     interfaces/InterfaceTypes.hpp
     interfaces/LinearOperator.hpp
     interfaces/MatrixBase.hpp
     interfaces/VectorBase.hpp
     solvers/KrylovSolver.hpp
     solvers/CGsolver.hpp
     solvers/BiCGSTABsolver.hpp
     solvers/PreconditionerBase.hpp
     solvers/PreconditionerIdentity.hpp
     utilities/BlockVectorView.hpp
     utilities/BlockVectorWrapper.hpp
     utilities/BlockVector.hpp
     utilities/BlockOperatorView.hpp
     utilities/BlockOperatorWrapper.hpp
     utilities/BlockOperator.hpp
     utilities/LinearSolverParameters.hpp
     utilities/LAIHelperFunctions.hpp
     DofManager.hpp
     DofManagerHelpers.hpp )

#
# Specify all sources
#
set( linearAlgebra_sources
     interfaces/BlasLapackLA.cpp
     solvers/CGsolver.cpp
     solvers/BiCGSTABsolver.cpp
     utilities/LAIHelperFunctions.cpp
     DofManager.cpp )

if( BUILD_OBJ_LIBS )
  set( dependencyList common blas lapack raja managers)
else()
  set( dependencyList common blas lapack raja )
endif()


if( ENABLE_MPI )
    list( APPEND dependencyList mpi )
endif()

if( ENABLE_TRILINOS )

    list( APPEND linearAlgebra_headers
          interfaces/trilinos/EpetraMatrix.hpp
          interfaces/trilinos/EpetraVector.hpp
          interfaces/trilinos/EpetraUtils.hpp
          interfaces/trilinos/TrilinosSolver.hpp
          interfaces/trilinos/TrilinosInterface.hpp )

    list( APPEND linearAlgebra_sources
          interfaces/trilinos/EpetraMatrix.cpp
          interfaces/trilinos/EpetraVector.cpp
          interfaces/trilinos/TrilinosSolver.cpp )

    list( APPEND dependencyList trilinos )

endif()

if( ENABLE_HYPRE )

    list( APPEND linearAlgebra_headers
          interfaces/hypre/HypreMatrix.hpp
          interfaces/hypre/HypreVector.hpp
          interfaces/hypre/HypreSolver.hpp
          interfaces/hypre/HypreInterface.hpp
          interfaces/hypre/HypreUtils.hpp )

    list( APPEND linearAlgebra_sources
          interfaces/hypre/HypreMatrix.cpp
          interfaces/hypre/HypreVector.cpp
          interfaces/hypre/HypreSolver.cpp )

    list( APPEND dependencyList hypre )

endif()

if( ENABLE_PETSC )

    list( APPEND linearAlgebra_headers
          interfaces/petsc/PetscVector.hpp
          interfaces/petsc/PetscMatrix.hpp
          interfaces/petsc/PetscSolver.hpp
          interfaces/petsc/PetscInterface.hpp
          interfaces/petsc/PetscUtils.hpp )

    list( APPEND linearAlgebra_sources
          interfaces/petsc/PetscVector.cpp
          interfaces/petsc/PetscMatrix.cpp
          interfaces/petsc/PetscSolver.cpp
          interfaces/petsc/PetscInterface.cpp )

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
                 HEADERS               ${linearAlgebra_headers}
                 DEPENDS_ON            ${dependencyList}
               )

target_include_directories( linearAlgebra PUBLIC ${CMAKE_CURRENT_LIST_DIR} )

geosx_add_code_checks( PREFIX linearAlgebra )

add_subdirectory( unitTests )

message( "Leaving /src/coreComponents/linearAlgebra/CMakeLists.txt")
