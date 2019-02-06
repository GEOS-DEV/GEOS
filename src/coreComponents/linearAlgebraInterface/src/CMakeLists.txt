#
# Specify all headers
#
set( linearAlgebraInterface_headers
     BlockMatrixView.hpp
     BlockVectorView.hpp
     DofManager.hpp
     EpetraMatrix.hpp
     EpetraVector.hpp
     HypreInterface.hpp 
     InterfaceTypes.hpp
     LinearSolverParameters.hpp
     NativeCG.hpp
     NativeBICGSTAB.hpp
     TrilinosSolver.hpp
     TrilinosInterface.hpp
     HypreInterface.hpp 
     LinearSolverParameters.hpp
     LapackInterface.hpp
     LapackMatrix.hpp)

#
# Specify all sources
#
set( linearAlgebraInterface_sources
     BlockVectorView.cpp
     BlockMatrixView.cpp
     DofManager.cpp
     EpetraMatrix.cpp
     EpetraVector.cpp
     TrilinosSolver.cpp
     BlockVectorView.cpp
     BlockMatrixView.cpp
     LapackMatrix.cpp)

set( dependencyList common trilinos hypre lapack_suite mpi)

if( CONDUIT_FOUND )
  set( dependencyList ${dependencyList} conduit )
endif()

blt_add_library( NAME                  linearAlgebraInterface
                 SOURCES               ${linearAlgebraInterface_sources}
                 HEADERS               ${linearAlgebraInterface_headers}
                 DEPENDS_ON            ${dependencyList}
               )
               
target_include_directories( linearAlgebraInterface PUBLIC ${CMAKE_CURRENT_LIST_DIR})

geosx_add_code_checks(PREFIX linearAlgebraInterface )
