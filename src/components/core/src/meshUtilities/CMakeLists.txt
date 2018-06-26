
#
# Specify all headers
#
set(meshUtilities_headers
    ComputationalGeometry.hpp
    MeshManager.hpp
    MeshGeneratorBase.hpp
    InternalMeshGenerator.hpp
    MeshUtilities.hpp
    SimpleGeometricObjects/GeometricObjectManager.hpp
    SimpleGeometricObjects/SimpleGeometricObjectBase.hpp
    SimpleGeometricObjects/Box.hpp
    StructuredGridUtilities.hpp
   )

#
# Specify all sources
#
set(meshUtilities_sources
    ComputationalGeometry.cpp
    MeshManager.cpp
    MeshGeneratorBase.cpp
    InternalMeshGenerator.cpp
    MeshUtilities.cpp
    SimpleGeometricObjects/GeometricObjectManager.cpp
    SimpleGeometricObjects/SimpleGeometricObjectBase.cpp
    SimpleGeometricObjects/Box.cpp
   )
    
blt_add_library( NAME                  meshUtilities
                 SOURCES               ${meshUtilities_sources}
                 HEADERS               ${meshUtilities_headers}
                 DEPENDS_ON            fileIO common trilinos
               )
               
target_include_directories( meshUtilities PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src)
geosx_add_code_checks(PREFIX meshUtilities )