#
# Specify all headers
#
set( linearAlgebraInterface_headers
     BlockMatrixView.hpp
     EpetraSparseMatrix.hpp
     EpetraVector.hpp
     TrilinosSolver.hpp
     TrilinosInterface.hpp )

#
# Specify all sources
#
set( linearAlgebraInterface_sources
     EpetraSparseMatrix.cpp
     EpetraVector.cpp
     TrilinosSolver.cpp
     TrilinosInterface.cpp)

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
