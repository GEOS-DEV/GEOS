
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    MeshUtilities/MeshManager.hpp
    MeshUtilities/MeshGeneratorBase.hpp
    MeshUtilities/InternalMeshGenerator.hpp
    MeshUtilities/MeshUtilities.hpp
    MeshUtilities/SimpleGeometricObjects/GeometricObjectManager.hpp
    MeshUtilities/SimpleGeometricObjects/SimpleGeometricObjectBase.hpp
    MeshUtilities/SimpleGeometricObjects/Box.hpp
    MeshUtilities/StructuredGridUtilities.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
    ${geosx_core_sources}
    MeshUtilities/MeshManager.cpp
    MeshUtilities/MeshGeneratorBase.cpp
    MeshUtilities/InternalMeshGenerator.cpp
    MeshUtilities/MeshUtilities.cpp
    MeshUtilities/SimpleGeometricObjects/GeometricObjectManager.cpp
    MeshUtilities/SimpleGeometricObjects/SimpleGeometricObjectBase.cpp
    MeshUtilities/SimpleGeometricObjects/Box.cpp
    PARENT_SCOPE )
