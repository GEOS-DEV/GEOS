
#
# Specify all headers
#
set(meshUtilities_headers
    ComputationalGeometry.hpp
    MeshManager.hpp
    MeshGeneratorBase.hpp
    InternalMeshGenerator.hpp
    InternalWellGenerator.hpp
    PerforationData.hpp
    Perforation.hpp
    MeshUtilities.hpp
    SimpleGeometricObjects/GeometricObjectManager.hpp
    SimpleGeometricObjects/SimpleGeometricObjectBase.hpp
    SimpleGeometricObjects/Box.hpp
    SimpleGeometricObjects/ThickPlane.hpp
    SimpleGeometricObjects/BoundedPlane.hpp
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
    InternalWellGenerator.cpp
    PerforationData.cpp
    Perforation.cpp
    MeshUtilities.cpp
    SimpleGeometricObjects/GeometricObjectManager.cpp
    SimpleGeometricObjects/SimpleGeometricObjectBase.cpp
    SimpleGeometricObjects/Box.cpp
    SimpleGeometricObjects/Cylinder.cpp
    SimpleGeometricObjects/ThickPlane.cpp
    SimpleGeometricObjects/BoundedPlane.cpp
   )

if( BUILD_OBJ_LIBS)
  set( dependencyList common dataRepository fileIO mesh)
else()
  set( dependencyList 
       common
       mesh
       mpiCommunications
       managers
       codingUtilities
       dataRepository
      )
endif()

if( ENABLE_PAMELA )
    message(STATUS "Adding PAMELAMeshGenerator sources and headers")
    set( meshUtilities_headers ${meshUtilities_headers} PAMELAMeshGenerator.hpp )
    set( meshUtilities_sources ${meshUtilities_sources} PAMELAMeshGenerator.cpp)
    set( dependencyList ${dependencyList} PAMELA )
endif()

if( ENABLE_OPENMP )
    set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  meshUtilities
                 SOURCES               ${meshUtilities_sources}
                 HEADERS               ${meshUtilities_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )
               
target_include_directories( meshUtilities PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

add_subdirectory( unitTests )

geosx_add_code_checks(PREFIX meshUtilities )
