#
# Specify all headers
#
set( codingUtilities_headers
     EnumStrings.hpp
     Parsing.hpp
     SFINAE_Macros.hpp
     StringUtilities.hpp
     UnitTestUtilities.hpp
     Utilities.hpp
     traits.hpp
   )

#
# Specify all sources
#
set( codingUtilities_sources
     StringUtilities.cpp )

set( dependencyList common fast_float )

if( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME       codingUtilities
                 SOURCES    ${codingUtilities_sources}
                 HEADERS    ${codingUtilities_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOSX_BUILD_OBJ_LIBS}
               )

target_include_directories( codingUtilities PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

geosx_add_code_checks(PREFIX codingUtilities )

add_subdirectory(tests)
