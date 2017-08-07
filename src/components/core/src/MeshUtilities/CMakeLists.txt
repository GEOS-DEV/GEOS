
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    MeshUtilities/MeshGenerator.hpp
    MeshUtilities/MeshUtilities.hpp
    MeshUtilities/SimpleGeometricObjects/SimpleGeometricObjectBase.hpp
    MeshUtilities/SimpleGeometricObjects/Box.hpp
    MeshUtilities/StructuredGridUtilities.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
    ${geosx_core_sources}
    MeshUtilities/MeshGenerator.cpp
    MeshUtilities/MeshUtilities.cpp
    MeshUtilities/SimpleGeometricObjects/SimpleGeometricObjectBase.cpp
    MeshUtilities/SimpleGeometricObjects/Box.cpp
    PARENT_SCOPE )
