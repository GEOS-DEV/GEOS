#
# Specify all headers
#
set( linearAlgebraInterface_headers
     NativeCG.hpp
     NativeBICGSTAB.hpp
     BlockMatrixView.hpp
     BlockVectorView.hpp
     EpetraMatrix.hpp
     EpetraVector.hpp
     InterfaceTypes.hpp
     TrilinosSolver.hpp
     TrilinosInterface.hpp
     HypreInterface.hpp 
     LinearSolverParameters.hpp
     LapackSuiteInterface.hpp
     BlasMatrix.hpp)

#
# Specify all sources
#
set( linearAlgebraInterface_sources
     EpetraMatrix.cpp
     EpetraVector.cpp
     TrilinosSolver.cpp
     BlockVectorView.cpp
     BlockMatrixView.cpp
     BlasMatrix.cpp)

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
