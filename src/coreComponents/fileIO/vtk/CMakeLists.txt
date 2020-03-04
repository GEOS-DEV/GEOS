#
# Specify all headers
#
set( vtk_interface_headers
     DummyVTKDataArray.hpp
     VerticesVTKDataArray.hpp )

#
# Specify all sources
#
set( vtk_interface_sources
     DummyVTKDataArray.cpp
     VerticesVTKDataArray.cpp )

set( dependencyList vtk )

blt_add_library( NAME                  vtk_interface
                 SOURCES               ${vtk_interface_sources}
                 HEADERS               ${vtk_interface_headers}
                 DEPENDS_ON            ${dependencyList}
               )
               
target_include_directories( vtk_interface PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX vtk_interface )
