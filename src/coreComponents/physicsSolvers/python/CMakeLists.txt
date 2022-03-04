set( pysolver_headers PySolver.hpp
  PySolverType.hpp
		    )

set( pysolver_sources PySolver.cpp
  )


blt_add_library( NAME             pysolver
                 SOURCES          ${pysolver_sources}
                 HEADERS          ${pysolver_headers}
                 DEPENDS_ON       ${dependencyList} ${extraComponentsList} ${externalComponentList} Python3::Python physicsSolvers
                 OBJECT           ${GEOSX_BUILD_OBJ_LIBS}
                 )

target_include_directories( pysolver PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

geosx_add_code_checks( PREFIX pysolver )
