#
# Specify all headers
#
set( linearAlgebraInterface_headers
     BlasLapackFunctions.h
     BlasLapackLA.hpp
     BlockMatrixView.hpp
     BlockVectorView.hpp
     DofManager.hpp
     InterfaceTypes.hpp
     LinearSolverParameters.hpp
     NativeCG.hpp
     NativeBICGSTAB.hpp
     TrilinosSolver.hpp
     TrilinosInterface.hpp
     PetscVector.hpp
     PetscSparseMatrix.hpp
     PetscSolver.hpp
     PetscInterface.hpp ) 

#
# Specify all sources
#
set( linearAlgebraInterface_sources
     BlasLapackLA.cpp
     BlockVectorView.cpp
     BlockMatrixView.cpp
     DofManager.cpp )

set( dependencyList common blas lapack raja hdf5 pugixml mpi )

if( ENABLE_TRILINOS )

  list( APPEND linearAlgebraInterface_headers
        EpetraMatrix.hpp
        EpetraVector.hpp
        TrilinosSolver.hpp
        TrilinosInterface.hpp )

  list( APPEND linearAlgebraInterface_sources
        EpetraMatrix.cpp
        EpetraVector.cpp
        TrilinosSolver.cpp )

endif()

if( ENABLE_HYPRE )

  list( APPEND linearAlgebraInterface_headers
        HypreInterface.hpp )

endif()

if( ENABLE_PETSC )

  list( APPEND linearAlgebraInterface_headers
        PetscVector.hpp
        PetscSparseMatrix.hpp
        PetscSolver.hpp
        PetscInterface.hpp )
      
  list( APPEND linearAlgebraInterface_sources 
        PetscVector.cpp
        PetscSparseMatrix.cpp
        PetscSolver.cpp )

  list( APPEND dependencyList petsc )

endif()

if (CONDUIT_FOUND)
  set (dependencyList ${dependencyList} conduit)
endif()

if( ENABLE_OPENMP )
    set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  linearAlgebraInterface
                 SOURCES               ${linearAlgebraInterface_sources}
                 HEADERS               ${linearAlgebraInterface_headers}
                 DEPENDS_ON            ${dependencyList}
                 #OBJECT                TRUE
               )
               
target_include_directories( linearAlgebraInterface PUBLIC ${CMAKE_CURRENT_LIST_DIR})

geosx_add_code_checks( PREFIX linearAlgebraInterface )
