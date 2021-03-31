#
# Specify all headers
#
set(mainInterface_headers
    NumericalMethodsManager.hpp
   )


#
# Specify all sources
#

set(mainInterface_sources
    NumericalMethodsManager.cpp
   )

blt_add_library( NAME                  discretizationMethods
                 SOURCES               ${mainInterface_sources}
                 HEADERS               ${mainInterface_headers}
                 DEPENDS_ON            finiteElement finiteVolume virtualElement
                 OBJECT                ${buildAsObj}
               )

target_include_directories( discretizationMethods PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks( PREFIX discretizationMethods )
