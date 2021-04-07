#
# Specify all sources
#
set( codingUtilities_sources
     StringUtilities.cpp
   )

set( dependencyList common )

blt_add_library( NAME                  codingUtilities
                 SOURCES               ${codingUtilities_sources}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )

target_include_directories( codingUtilities PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX codingUtilities )

add_subdirectory(tests)
