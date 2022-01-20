set( pysolver_headers PySolver.hpp
  PySolverType.hpp
		    )

set( pysolver_sources PySolver.cpp
  )

set( dependencyList ${extraComponentsLinkList} ${externalComponentsLinkList} Python3::Python )

blt_add_library( NAME             pysolver
                 SOURCES          ${pysolver_sources}
                 HEADERS          ${pysolver_headers}
                 DEPENDS_ON       physicsSolvers ${dependencyList}
                 SHARED TRUE
                 CLEAR_PREFIX TRUE
                 )

target_include_directories( pysolver PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

geosx_add_code_checks( PREFIX pysolver )
