#
# Specify all headers
#
set( denseLinearAlgebra_headers
     interfaces/blaslapack/BlasLapackFunctions.h
     interfaces/blaslapack/BlasLapackLA.hpp
   )
#
# Specify all sources
#
set( denseLinearAlgebra_sources
     interfaces/blaslapack/BlasLapackLA.cpp
   )

set( dependencyList codingUtilities blas lapack )

if ( ENABLE_CUDA )
    list( APPEND dependencyList cuda )
endif()

blt_add_library( NAME       denseLinearAlgebra
                 SOURCES    ${denseLinearAlgebra_sources}
                 HEADERS    ${denseLinearAlgebra_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOSX_BUILD_OBJ_LIBS}
               )

target_include_directories( denseLinearAlgebra PUBLIC ${CMAKE_CURRENT_LIST_DIR} )

add_subdirectory( unitTests )

geosx_add_code_checks( PREFIX denseLinearAlgebra )
