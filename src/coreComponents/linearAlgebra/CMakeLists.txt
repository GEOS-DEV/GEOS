message( "Entering /src/coreComponents/linearAlgebra/CMakeLists.txt")

#
# Specify all headers
#
set( linearAlgebra_headers
     interfaces/BlasLapackFunctions.h
     interfaces/BlasLapackLA.hpp
     interfaces/InterfaceTypes.hpp
     solvers/CGsolver.hpp
     solvers/BiCGSTABsolver.hpp
     utilities/BlockMatrixView.hpp
     utilities/BlockVectorView.hpp
     utilities/LinearSolverParameters.hpp
     DofManager.hpp
     DofManagerHelpers.hpp
     utilities/LAIHelperFunctions.hpp)

#
# Specify all sources
#
set( linearAlgebra_sources
     interfaces/BlasLapackLA.cpp
     solvers/CGsolver.cpp
     solvers/BiCGSTABsolver.cpp
     utilities/BlockVectorView.cpp
     utilities/BlockMatrixView.cpp
     DofManager.cpp
     utilities/LAIHelperFunctions.cpp)

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
                 interfaces/EpetraMatrix.hpp
                 interfaces/EpetraVector.hpp
                 interfaces/TrilinosSolver.hpp
                 interfaces/TrilinosInterface.hpp )

    list( APPEND linearAlgebra_sources
                 interfaces/EpetraMatrix.cpp
                 interfaces/EpetraVector.cpp
                 interfaces/TrilinosSolver.cpp )

    list( APPEND dependencyList trilinos )

endif()

if( ENABLE_HYPRE )

    list( APPEND linearAlgebra_headers
                 interfaces/HypreInterface.hpp )

    list( APPEND dependencyList hypre )

endif()

if( ENABLE_PETSC )

    list( APPEND linearAlgebra_headers
                 interfaces/PetscVector.hpp
                 interfaces/PetscSparseMatrix.hpp
                 interfaces/PetscSolver.hpp
                 interfaces/PetscInterface.hpp )

    list( APPEND linearAlgebra_sources
                 interfaces/PetscVector.cpp
                 interfaces/PetscSparseMatrix.cpp
                 interfaces/PetscSolver.cpp
                 interfaces/PetscInterface.cpp)

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
