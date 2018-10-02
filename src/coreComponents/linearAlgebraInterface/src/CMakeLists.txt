#
# Specify all headers
#
set( linearAlgebraInterface_headers
     solver_CG.hpp
     solver_BiCGSTAB.hpp
     BlockMatrixView.hpp
     BlockVectorView.hpp
     TrilinosSparseMatrix.hpp
     TrilinosVector.hpp
     InterfaceTypes.hpp
     TrilinosSolver.hpp
     TrilinosInterface.hpp
     HypreInterface.hpp )

#
# Specify all sources
#
set( linearAlgebraInterface_sources
     TrilinosSparseMatrix.cpp
     TrilinosVector.cpp
     TrilinosSolver.cpp)

set( dependencyList common trilinos hypre mpi)

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
