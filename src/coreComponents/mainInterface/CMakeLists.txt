#
# Specify all headers
#
set( mainInterface_headers
     GeosxState.hpp
     ProblemManager.hpp
     initialization.hpp
     DumpToJson.hpp
   )

#
# Specify all sources
#
set( mainInterface_sources
     GeosxState.cpp
     ProblemManager.cpp
     initialization.cpp
   )

set( dependencyList physicsSolvers discretizationMethods fieldSpecification linearAlgebra events fileIO optionparser )


include(FetchContent)

FetchContent_Declare( json
                      GIT_REPOSITORY https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent
                      GIT_TAG v3.10.5
                      )
FetchContent_MakeAvailable( json )

set( dependencyList nlohmann_json::nlohmann_json ${dependencyList} )


if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME       mainInterface
                 SOURCES    ${mainInterface_sources}
                 HEADERS    ${mainInterface_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOSX_BUILD_OBJ_LIBS}
               )

add_dependencies( mainInterface generate_version )

target_include_directories( mainInterface PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

geosx_add_code_checks( PREFIX mainInterface )
