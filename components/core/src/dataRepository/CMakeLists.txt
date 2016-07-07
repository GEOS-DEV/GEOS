
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    dataRepository/DataTypes.hpp
    dataRepository/intrinsic/DataObject.hpp
    dataRepository/intrinsic/DataObjectBase.hpp
    dataRepository/intrinsic/DataObjectManager.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    dataRepository/intrinsic/DataObject.cpp
    dataRepository/intrinsic/DataObjectBase.cpp
    dataRepository/intrinsic/DataObjectManager.cpp
    PARENT_SCOPE )
