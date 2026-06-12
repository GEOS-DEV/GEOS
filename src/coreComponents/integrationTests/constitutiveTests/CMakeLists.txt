set( tplDependencyList ${parallelDeps} gtest )

set( dependencyList mainInterface  )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )
