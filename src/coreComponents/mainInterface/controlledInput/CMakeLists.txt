include( FetchContent )
FetchContent_Declare( yaml-cpp
                      GIT_REPOSITORY https://github.com/jbeder/yaml-cpp
                      GIT_TAG 0.8.0
)
FetchContent_MakeAvailable( yaml-cpp )

set( controlledInput_headers
     ControlledInput.hpp
     Events.hpp
     Mesh.hpp
     Misc.hpp
     Outputs.hpp
     Solvers.hpp
)

set( controlledInput_sources
     ControlledInput.cpp
     Events.cpp
     Mesh.cpp
     Misc.cpp
     Outputs.cpp
     Solvers.cpp
)

set( dependencyList dataRepository yaml-cpp )

blt_add_library( NAME       controlledInput
                 SOURCES    ${controlledInput_sources}
                 HEADERS    ${controlledInput_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOSX_BUILD_OBJ_LIBS}
)
