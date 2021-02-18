
#
# Specify all headers
#
set( codingUtilities_headers
     traits.hpp
     EnumStrings.hpp
     StringUtilities.hpp
     UnitTestUtilities.hpp
     Utilities.hpp
     static_if.hpp
     SFINAE_Macros.hpp
   )

#
# Specify all sources
#
set( codingUtilities_sources
     StringUtilities.cpp
   )

set( dependencyList common )

blt_add_library( NAME                  codingUtilities
                 SOURCES               ${codingUtilities_sources}
                 HEADERS               ${codingUtilities_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )

target_include_directories( codingUtilities PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX codingUtilities )

add_subdirectory(tests)
