
#add_subdirectory(src)
#add_subdirectory(unitTests)

#
# Specify all headers
#
set( linearAlgebraInterface_headers
     src/EpetraSparseMatrix.hpp
     src/EpetraVector.hpp
     src/TrilinosSolver.hpp
     src/TrilinosInterface.hpp )

#
# Specify all sources
#
set( linearAlgebraInterface_sources
     src/EpetraSparseMatrix.cpp
     src/EpetraVector.cpp
     src/TrilinosSolver.cpp
     src/TrilinosInterface.cpp)

set( dependencyList common trilinos hypre mpi)

blt_add_library( NAME                  linearAlgebraInterface
                 SOURCES               ${linearAlgebraInterface_sources}
                 HEADERS               ${linearAlgebraInterface_headers}
                 DEPENDS_ON            ${dependencyList}
               )
               
target_include_directories( linearAlgebraInterface PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src)

geosx_add_code_checks(PREFIX linearAlgebraInterface )
