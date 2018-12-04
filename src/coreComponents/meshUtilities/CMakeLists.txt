
#
# Specify all headers
#
set(meshUtilities_headers
    ComputationalGeometry.hpp
    MeshManager.hpp
    MeshGeneratorBase.hpp
    InternalMeshGenerator.hpp
    VTMMeshGenerator.hpp
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
    VTMMeshGenerator.cpp
    MeshUtilities.cpp
    SimpleGeometricObjects/GeometricObjectManager.cpp
    SimpleGeometricObjects/SimpleGeometricObjectBase.cpp
    SimpleGeometricObjects/Box.cpp
   )


set( dependencyList fileIO common trilinos )

set(ENABLE_PAMELA ON)
if( ENABLE_PAMELA )
    message(STATUS "Adding PAMELAMeshGenerator sources and headers")
    set( meshUtilities_headers ${meshUtilities_headers} PAMELAMeshGenerator.hpp )
    set( meshUtilities_sources ${meshUtilities_sources} PAMELAMeshGenerator.cpp)
    set( dependencyList ${dependencyList} PAMELA )
endif()

    
blt_add_library( NAME                  meshUtilities
                 SOURCES               ${meshUtilities_sources}
                 HEADERS               ${meshUtilities_headers}
                 DEPENDS_ON            ${dependencyList}
               )
               
target_include_directories( meshUtilities PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX meshUtilities )
