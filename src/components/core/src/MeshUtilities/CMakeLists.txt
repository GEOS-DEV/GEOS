
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    MeshUtilities/MeshGenerator.hpp
    MeshUtilities/SimpleGeometricObjects/SimpleGeometricObjectBase.hpp
    MeshUtilities/StructuredGridUtilities.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
    ${geosx_core_sources}
    MeshUtilities/MeshGenerator.cpp
    MeshUtilities/SimpleGeometricObjects/SimpleGeometricObjectBase.cpp
    PARENT_SCOPE )
