
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    dataRepository/DataObject.hpp
    dataRepository/DataObjectBase.hpp
    dataRepository/DataObjectManager.hpp
    dataRepository/DataTypes.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    dataRepository/DataObject.cpp
    dataRepository/DataObjectBase.cpp
    dataRepository/DataObjectManager.cpp
    PARENT_SCOPE )
