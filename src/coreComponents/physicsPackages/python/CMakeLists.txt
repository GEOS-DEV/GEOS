
if( ENABLE_PYGEOSX )
  set( physicsPackages_headers
       ${physicsPackages_headers}
       python/PySolver.hpp
       python/PySolverType.hpp
       PARENT_SCOPE )
  set( physicsPackages_sources
       ${physicsPackages_sources}
       python/PySolver.cpp
       PARENT_SCOPE )

  set( dependencyList ${dependencyList} Python3::Python pylvarray PARENT_SCOPE)
endif()

