include( FetchContent )
FetchContent_Declare( yaml-cpp
                      GIT_REPOSITORY https://github.com/jbeder/yaml-cpp
                      GIT_TAG 0.8.0
)
FetchContent_MakeAvailable( yaml-cpp )

set( controlledInput_headers
     BoundaryConditions.hpp
     ControlledInput.hpp
     Events.hpp
     Functions.hpp
     Geometry.hpp
     Mesh.hpp
     NumericalStrategies.hpp
     Misc.hpp
     Outputs.hpp
     Simulation.hpp
     Solvers.hpp
)

set( controlledInput_sources
     BoundaryConditions.cpp
     ControlledInput.cpp
     Events.cpp
     Functions.cpp
     Geometry.cpp
     Mesh.cpp
     NumericalStrategies.cpp
     Misc.cpp
     Outputs.cpp
     Simulation.cpp
     Solvers.cpp
)

set( dependencyList dataRepository yaml-cpp )

blt_add_library( NAME       controlledInput
                 SOURCES    ${controlledInput_sources}
                 HEADERS    ${controlledInput_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOSX_BUILD_OBJ_LIBS}
)
